library(shiny)
library(leaflet)
library(dplyr)

ui <- navbarPage('Alpha diversity', id = 'nav', 
                 tabPanel('Map', 
                          div(class = 'outer',
                              tags$head(
                                includeCSS('/Users/taylorminich/Documents/grad/fall 2018/genomics - 750/final/style.css'),
                                includeScript('/Users/taylorminich/Documents/grad/fall 2018/genomics - 750/final/map.js')
                              ),
                              leafletOutput('map', width = '100%', height = '100%'),
                              absolutePanel(id = 'controls', class = 'panel panel-default', fixed = TRUE,
                                            draggable = TRUE, top = 60, left = 'auto', right = 20, bottom = 'auto',
                                            width = 700, height = 'auto',
                                            h2('Site explorer'),
                                            checkboxInput('group', 'Group by site', TRUE),
                                            checkboxInput('regression', 'Add regression line', TRUE),
                                            checkboxGroupInput('site', 'site', c('Bartlett Experimental Forest' = 'bart',
                                                                                 'Central Plains Experimental Range' = 'cper',
                                                                                 'Disney Wilderness Preserve' = 'dsny',
                                                                                 'Harvard Forest' = 'harv',
                                                                                 'Jones Ecological Research Center' = 'jerc',
                                                                                 'Ordway-Swisher Biological Station' = 'osbs',
                                                                                 'Smithsonian Conservation Biology Institute' = 'scbi',
                                                                                 'North Sterling' = 'ster',
                                                                                 'Talladega National Forest' = 'tall',
                                                                                 'Woodworth' = 'wood'), 
                                                               selected = c('bart','cper','dsny','harv','jerc','osbs','scbi','ster','tall','wood')),
                                            plotOutput('plot', height = 300)))))
                                            
server <- function(input, output, session) {
  samples <- read.csv('/Users/taylorminich/Documents/grad/fall 2018/genomics - 750/final/samples.csv')
  sites <- read.csv('/Users/taylorminich/Documents/grad/fall 2018/genomics - 750/final/sites.csv')
  alpha <- read.csv('/Users/taylorminich/Documents/grad/fall 2018/genomics - 750/final/alpha.csv')
  samples <- merge(samples, alpha, by.x='id', by.y='sample')
  
  df <- reactive({
   samples %>%
      filter(site %in% input$site)

  })
  output$plot <- renderPlot({
    if(input$group == TRUE & input$regression == TRUE){
      ggplot(df(), aes(x = structure, y = func, color = site)) + 
        geom_point() + 
        geom_smooth(method = lm) +
        xlim(5, 6) + xlab('Structural diversity') +
        ylim(2.8, 2.95) + ylab('Functional diversity') 
    } else if(input$group == TRUE & input$regression == FALSE){
      ggplot(df(), aes(x = structure, y = func, color = site)) + 
        geom_point() + 
        xlim(5, 6) + xlab('Structural diversity') +
        ylim(2.8, 2.95) + ylab('Functional diversity') 
    } else if(input$group == FALSE & input$regression == TRUE){
      ggplot(df(), aes(x = structure, y = func)) + 
        geom_point() + 
        geom_smooth(method = lm) +
        xlim(5, 6) + xlab('Structural diversity') +
        ylim(2.8, 2.95) + ylab('Functional diversity') 
    } else {
      ggplot(df(), aes(x = structure, y = func)) + 
        geom_point() +
        xlim(5, 6) + xlab('Structural diversity') +
        ylim(2.8, 2.95) + ylab('Functional diversity') 
      }
  })
  output$map <- renderLeaflet({
    leaflet() %>% setView(lng = -60, lat = 38, zoom = 4) %>% addTiles()%>%
      addCircleMarkers(data = df(),~lon, ~lat, label = ~as.character(id), radius = 6, color = 'red', fillOpacity = 0.5) %>%
      addMarkers(data = sites, ~lon, ~lat, label = ~as.character(name))
  })

}

shinyApp(ui = ui, server = server)




##### Load in necessary packages #####

library(shiny)
library(ggplot2)
library(dplyr)
library(plotly)
library(shinythemes)
library(ggpubr)
library(shinydashboard)


ui <- dashboardPage(skin = 'black',
                    dashboardHeader(title = 'Microglia Development RNASeq',
                                    titleWidth = 400,
                                    dropdownMenu(type = "notifications",
                                                 notificationItem(
                                                   text = "WARNING LOW EXPRESSION",
                                                   icon("users")))),
                    dashboardSidebar(
                      sidebarMenu(
                        menuItem('Plots', tabName = 'Plots', icon = icon('bar-chart-o'), startExpanded = FALSE,
                                 menuSubItem('Bar Plot', tabName = 'BarPlot'),
                                 menuSubItem('Violin Plot', tabName = 'ViolinPlot'),
                                 menuSubItem('Line Plot', tabName = 'LinePlot'),
                                 menuSubItem('Dot Plot', tabName = 'DotPlot')),
                        menuItem('Info', tabName = "SummaryDataTable", icon = icon('info')),
                        menuItem('Stats', tabName = 'ANOVATables', icon = icon('th'))
                      )),
                    
                    dashboardBody(
                      tabItems(
                        tabItem(tabName = 'BarPlot', h2('Bar Plot goes here')),
                        tabItem(tabName = 'ViolinPlot', h2('Violin Plot goes here')),
                        tabItem('DotPlot', h2('Dot plot goes here')),
                        tabItem('LinePlot', h2('Line plot goes here')),
                        tabItem('SummaryDataTable', h2('Summary Table goes here')),
                        tabItem(tabName = 'ANOVATables', h2('ANOVA table goes here'))
                    )))




server <- function (input, output) {
  
}

shinyApp(ui = ui, server = server)

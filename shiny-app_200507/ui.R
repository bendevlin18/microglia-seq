

######################################
###                                ###
###       ui.R file for            ###
###       microglia-seq            ###
###        application             ###
###                                ###
######################################


ui <- dashboardPage(skin = 'black',
                  
                  dashboardHeader(title = 'Microglia Development RNASeq', titleWidth = 375),
                  
                  dashboardSidebar(
                    selectizeInput('entered_genes', multiple = TRUE, choices = NULL, label = 'Search for one or multiple genes',
                                   options = list(maxItems = 5, placeholder = 'Example: Cx3cr1')),
                    sidebarMenu(
                      menuItem('Plots', tabName = 'Plots', icon = icon('bar-chart-o'), startExpanded = FALSE,
                               menuSubItem('Bar Plot', tabName = 'BarPlot'),
                               menuSubItem('Violin Plot', tabName = 'ViolinPlot'),
                               menuSubItem('Line Plot', tabName = 'LinePlot'),
                               menuSubItem('Dot Plot', tabName = 'DotPlot')),
                      menuItem('Info', tabName = "SummaryDataTable", icon = icon('info')),
                      menuItem('Stats', tabName = 'ANOVATable', icon = icon('th'))
                    ),
                    
                    checkboxGroupInput('checked_groups', label = 'Choose two ages to compare (Two-Way Anova)', 
                                       choices = c('E18', 'P4', 'P14', 'P60', 'P60 + LPS')),
                    checkboxInput('twoway_anova', label = 'Compare Across All Ages (Two-Way Anova)'),
                    checkboxInput('p_vals', label = 'Turn P values ON (For Bar, Violin, and single gene Line plots ONLY)'),
                    checkboxInput('ind_points', label = 'Turn Individual Points ON')
                  ),
                  
                  dashboardBody(
                    tabItems(
                      tabItem(tabName = 'BarPlot', plotlyOutput('bar', width = 'auto', height = 'auto')),
                      tabItem(tabName = 'ViolinPlot', plotlyOutput('violin', width = 'auto', height = 'auto')),
                      tabItem('DotPlot', plotlyOutput('dotplot', width = 'auto', height = 'auto')),
                      tabItem('LinePlot', plotlyOutput('line', width = 'auto', height = 'auto')),
                      tabItem('SummaryDataTable', 
                              box(tableOutput('summary_table')),
                              box(htmlOutput('open_target_link')),
                              box(htmlOutput('open_target_link2')))
                        
                        ,
                      tabItem(tabName = 'ANOVATable', tableOutput('anova_table'))))
)





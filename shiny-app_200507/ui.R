

######################################
###                                ###
###       ui.R file for            ###
###       microglia-seq            ###
###        application             ###
###                                ###
######################################


ui <-  dashboardPage(skin = 'black', 
                  
                  dashboardHeader(title = 'Microglia Development RNASeq', titleWidth = 375),
                  
                  dashboardSidebar(
                    selectizeInput('entered_genes', multiple = TRUE, choices = NULL, label = 'Search for one or multiple genes',
                                   options = list(maxItems = 5, placeholder = 'Example: Cx3cr1')),
                    sidebarMenu(
                      menuItem('Landing Page/Welcome', tabName = 'Landing', icon = icon('brain')),
                      menuItem('Plots', tabName = 'Plots', icon = icon('bar-chart-o'), startExpanded = FALSE,
                               menuSubItem('Bar Plot', tabName = 'BarPlot', icon = icon('bar-chart-o')),
                               menuSubItem('Violin Plot', tabName = 'ViolinPlot', icon = icon('chart-area')),
                               menuSubItem('Line Plot', tabName = 'LinePlot', icon = icon('chart-line')),
                               menuSubItem('Dot Plot', tabName = 'DotPlot', icon = icon('circle'))),
                      menuItem('Info', tabName = "SummaryDataTable", icon = icon('info')),
                      menuItem('Stats', tabName = 'ANOVATable', icon = icon('divide')))),
                  
                  dashboardBody(
                    tabItems(
                      tabItem(tabName = 'Landing', htmlOutput('welcome_info')),
                      tabItem(tabName = 'BarPlot', 
                              tags$style(type = "text/css", "#bar {height: calc(100vh - 80px) !important;}"),
                              column(width = 8,
                                     box(width = NULL, plotlyOutput('bar', width = '100%', height = '90%'))),
                              column(width = 4,
                                     fluidRow(
                                       box(width = 4, checkboxInput('p_vals', label = 'Turn P values ON')),
                                       box(width = 4, checkboxInput('ind_points', label = 'Turn Individual Points ON')))))
                              
                              
                              
                              
                              
                              
                              
                              
                              
                              
                              
                              
                              ,
                      tabItem(tabName = 'ViolinPlot',
                              tags$style(type = "text/css", "#violin {height: calc(100vh - 80px) !important;}"),
                              plotlyOutput('violin', width = '80%', height = '90%'),
                              checkboxInput('p_vals', label = 'Turn P values ON'),
                              checkboxInput('ind_points', label = 'Turn Individual Points ON')),
                      tabItem('DotPlot', 
                              tags$style(type = "text/css", "#dotplot {height: calc(100vh - 80px) !important;}"),
                              plotlyOutput('dotplot', width = '80%', height = '90%')),
                      tabItem('LinePlot', 
                              tags$style(type = "text/css", "#line {height: calc(100vh - 80px) !important;}"),
                              plotlyOutput('line', width = '80%', height = '90%'),
                              checkboxInput('p_vals', label = 'Turn P values ON (For Single Gene Line plots ONLY)'),
                              checkboxInput('ind_points', label = 'Turn Individual Points ON')),
                      tabItem('SummaryDataTable', 
                              box(DTOutput('summary_table')),
                              box(htmlOutput('open_target_link'))),
                      tabItem(tabName = 'ANOVATable', 
                              box(DTOutput('anova_table'), width = '75%'),
                              box(checkboxInput('twoway_anova', label = 'Compare Across All Ages (Two-Way Anova)')),
                              box(checkboxGroupInput('checked_groups', label = 'Choose two ages to compare (Two-Way Anova)', 
                                                     choices = c('E18', 'P4', 'P14', 'P60', 'P60 + LPS'))))
                              )))






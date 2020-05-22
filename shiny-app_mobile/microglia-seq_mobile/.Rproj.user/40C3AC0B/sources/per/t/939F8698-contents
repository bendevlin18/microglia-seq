

######################################
###                                ###
###       ui.R file for            ###
###       microglia-seq            ###
###        application             ###
###                                ###
######################################


ui <-  dashboardPage(skin = 'black', 
                  
                  dashboardHeader(title = 'Microglia Development RNASeq', titleWidth = 350),
                  
                  dashboardSidebar(
                    width = 275,
                    sidebarMenu(
                      menuItem('Landing Page', tabName = 'Landing', icon = icon('brain')),
                      menuItem('MDI Plots', tabName = 'mdi_plots', icon = icon('bar-chart-o'), startExpanded = FALSE,
                               menuSubItem('MDI Across Ages', tabName = 'mdi', icon = icon('bar-chart-o')),
                               menuSubItem('Heatmap of Index Genes', tabName = 'mdi_heat', icon = icon('bars')),
                               menuSubItem('List of Index Genes', tabName = 'mdi_gene_list', icon = icon('clipboard-list'))),
                      menuItem('TPM Plots', tabName = 'Plots', icon = icon('bar-chart-o'), startExpanded = FALSE,
                               selectizeInput('entered_genes', multiple = TRUE, choices = NULL, label = 'Search for one or multiple genes',
                                              options = list(maxItems = 4, placeholder = 'Example: Cx3cr1')),
                               menuSubItem('Bar Plot', tabName = 'BarPlot', icon = icon('bar-chart-o')),      
                               menuSubItem('Violin Plot', tabName = 'ViolinPlot', icon = icon('chart-area')),
                               menuSubItem('Line Plot', tabName = 'LinePlot', icon = icon('chart-line')),
                               menuSubItem('Dot Plot', tabName = 'DotPlot', icon = icon('circle')),
                               menuSubItem('Data Table', tabName = "SummaryDataTable", icon = icon('info')),
                               menuSubItem('Stats', tabName = 'ANOVATable', icon = icon('divide')),
                               menuSubItem('What is TPM?', tabName = 'tpm_info', icon = icon('question'))),
                      menuItem('Website Information', tabName = 'web_info', icon = icon('laptop-code')))),
                  
                  dashboardBody(
                    tabItems(
                      tabItem(tabName = 'Landing',
                              class = 'text-center',
                              uiOutput('landing_text'),
                              tags$style("#landing_text{color: black; font-size: 15px;}"),
                              tags$iframe(height = '450px', width = '50%', src='https://www.youtube.com/embed/8UhJ4CvmNgg')),
                      ## https://community.rstudio.com/t/how-to-embed-videos-in-shiny/38937/2
                      
                      
                      tabItem(tabName = 'mdi', 
                              tags$style(type = "text/css", "#mdi_plot {height: calc(100vh - 80px) !important;}"),
                              column(width = 9,
                                     plotlyOutput('mdi_plot', width = '100%', height = '90%')),
                              column(width = 3,
                                     fluidRow(
                                       checkboxInput('p_vals_mdi', label = 'Turn Male/Female P values ON')))),
                      tabItem(tabName = 'mdi_heat', 
                              tags$style(type = "text/css", "#mdi_heatmap {height: calc(100vh - 80px) !important;}"),
                              column(width = 9,
                                     withSpinner(plotlyOutput('mdi_heatmap', width = '100%', height = '90%'))),
                              column(width = 3,
                                     fluidRow(
                                       numericInput('n_heatmap_genes', label = 'How many genes?', value = 5, min = 1, max = 50)))),
                      tabItem(tabName = 'mdi_gene_list',
                              DTOutput('mdi_genes'),
                              p(class = 'text-center', downloadButton('mdi_gene_download', 'MDI Gene Download'))),
                      
                      
                      tabItem(tabName = 'BarPlot', 
                              tags$style(type = "text/css", "#bar {height: calc(100vh - 80px) !important;}"),
                              column(width = 9,
                                     withSpinner(plotlyOutput('bar', width = '100%', height = '90%'))),
                              column(width = 3,
                                     fluidRow(
                                       checkboxInput('p_vals', label = 'Turn Male/Female P values ON'),
                                       checkboxInput('ind_points', label = 'Turn Individual Points ON')))),
                      tabItem(tabName = 'ViolinPlot',
                              tags$style(type = "text/css", "#violin {height: calc(100vh - 80px) !important;}"),
                              column(width = 10,
                                     withSpinner(plotlyOutput('violin', width = '100%', height = '90%'))),
                              column(width = 2,
                                     fluidRow(
                                       checkboxInput('p_vals_violin', label = 'Turn Male/Female P values ON'),
                                       checkboxInput('ind_points_violin', label = 'Turn Individual Points ON')))),
                      tabItem('DotPlot', 
                              tags$style(type = "text/css", "#dotplot {height: calc(100vh - 80px) !important;}"),
                              withSpinner(plotlyOutput('dotplot', width = '80%', height = '90%'))),
                      tabItem('LinePlot', 
                              tags$style(type = "text/css", "#line {height: calc(100vh - 80px) !important;}"),
                              column(width = 10,
                                     withSpinner(plotlyOutput('line', width = '100%', height = '90%'))),
                              column(width = 2,
                                     fluidRow(
                                       checkboxInput('p_vals_line', label = 'Turn Male/Female P values ON')))),
                      tabItem('SummaryDataTable', 
                              withSpinner(DTOutput('summary_table')),
                              p(class = 'text-center', downloadButton('summary_table_download', 'Download Summary Table')),
                              tags$head(tags$style("#open_target_link{color: black; font-size: 25px; font-style: bold;}"),
                                        tags$style("#open_target_link2{color: black; font-size: 25px; font-style: bold;}"),
                                        tags$style("#open_target_link3{color: black; font-size: 25px; font-style: bold;}"),
                                        tags$style("#open_target_link4{color: black; font-size: 25px; font-style: bold;}")),
                              class = 'text-center',
                              htmlOutput('open_target_link'),
                              htmlOutput('open_target_link2'),
                              htmlOutput('open_target_link3'),
                              htmlOutput('open_target_link4')),
                      tabItem(tabName = 'ANOVATable', 
                              withSpinner(DTOutput('anova_table')),
                              p(class = 'text-center', downloadButton('anova_download', 'Download Anova Table')),
                              checkboxGroupInput('checked_groups', label = 'Choose Ages to compare (Two-Way Anova)', 
                                                  choices = c('All Ages', 'E18', 'P4', 'P14', 'P60', 'P60 + LPS'),
                                                  selected = 'All Ages')),
                      tabItem(tabName = 'tpm_info',
                              class = 'text-center',
                              uiOutput('TPM_info'),
                              tags$style("#tpm_info{color: black; font-size: 15px;}"),
                              tags$iframe(height = '450px', width = '50%', src = "https://www.youtube.com/embed/TTUrtCY2k-w", 
                                          frameborder="0", allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture", 
                                          allowfullscreen = TRUE)),
                      tabItem(tabName = 'web_info',
                              class = 'text-center',
                              uiOutput('website_info_text'),
                              tags$style("#website_info_text{color: black;font-size: 15px}"))

                              )))






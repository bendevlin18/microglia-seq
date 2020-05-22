


################################
################################
#####                      #####
#####                      #####
#####      ui.r for        #####
#####   microglia-seq      #####
#####      mobile          #####
#####                      #####
#####                      #####
################################
################################


ui <-  f7Page(
  title = "Microglia-Seq Mobile-Friendly",
  init = f7Init(theme = "light", skin = "ios"),
  f7TabLayout(
    panels = tagList(
      f7Panel(side = "left", theme = "dark", effect = "reveal",
              selectizeInput('entered_genes', multiple = TRUE, choices = NULL, label = 'Search for a gene',
                             options = list(maxItems = 1, placeholder = 'Example: Cx3cr1'))),
      f7Panel(side = "right", theme = "dark", effect = "reveal",
              f7checkBox('p_vals', label = 'Turn Male/Female Significance Levels ON'),
              f7checkBox('ind_points', label = 'Turn Individual Points ON'))),
    navbar = f7Navbar(
      title = "Microglia-Seq Mobile", left_panel = TRUE, right_panel = TRUE),
    f7Tabs(
      animated = TRUE, id = 'tabs',
      f7Tab(tabName = "Landing", active = TRUE, icon = f7Icon('home_fill'),
            column(width = 12,
                   tags$style("#landing_text{color: black; font-size: 20px;}"),
                   htmlOutput('landing_text')),
            tags$iframe(height = '250px', width = '100%', src='https://www.youtube.com/embed/8UhJ4CvmNgg', frameborder="0", 
                        allow="accelerometer; autoplay; encrypted-media; gyroscope; picture-in-picture", allowfullscreen=NA)),
      f7Tab(tabName = "MDI Plot", active = FALSE, icon = f7Icon('layers'),
            column(width = 4,
                   withSpinner(plotlyOutput('mdi_plot'))),
            column(width = 4,
                   withSpinner(DTOutput('mdi_genes')))),
      f7Tab(tabName = "Bar Plot", active = FALSE, icon = f7Icon('graph_square', lib = 'ios'),
            column(width = 4,
                   withSpinner(plotlyOutput('bar'))),
            column(width = 3,
                   withSpinner(DTOutput('summary_table')))),
      f7Tab(tabName = "Line Plot", active = FALSE, icon = f7Icon('bars', lib = 'ios'),
            column(width = 4,
                   withSpinner(plotlyOutput('line'))),
            column(width = 12,
                   withSpinner(DTOutput('line_summary_table')))),
      f7Tab(tabName = "Dot Plot", active = FALSE, icon = f7Icon('graph_circle_fill'),
            column(width = 4,
                   withSpinner(plotlyOutput('dotplot'))),
            column(width = 6,
                   withSpinner(DTOutput('dot_summary_table')))),
      f7Tab(tabName = "Site Info", active = FALSE, icon = f7Icon('info_fill'),
            column(width = 4,
                   uiOutput('website_info_text'))))
))

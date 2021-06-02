library(shiny)
library(ggplot2)
library(dplyr)
library(shinydashboard)


##using shiny and R to emulate the interactive data interface seen here: https://www.brainrnaseq.org/

##hope to model some of the interactive bits of this code off of this app: https://yihanw.shinyapps.io/Recipe_Nutrition/?_ga=2.68778321.1739234079.1586842901-780141392.1586842901

##how to integrate shiny apps with squarespace websites


df <- as_tibble(read.csv('GSE99622_hanamsagar2017_cleaned_melted.csv'))

genes = unique(df['gene'])

std_err <- function(x) sd(x) / sqrt(length(x))

df_grouped <- df %>%
  group_by(gene, age, sex) %>%
  summarise(sem = std_err(expression), expression = mean(expression), n = n())



ui <- fluidPage(
                
  titlePanel('Microglia Development RNASeq'),
  
  sidebarLayout(
    sidebarPanel(
      selectizeInput('entered_gene', label = 'Search for a gene', choices = genes,
                   options = list(
                     placeholder = 'Example: Csf1r',
                     onInitialize = I('function() { this.setValue(""); }')))),
    mainPanel(
      tabsetPanel(
        tabPanel('Bar Plot', plotOutput('bar')),
        tabPanel('Violin Plot', plotOutput('violin')),
        tabPanel('Summary Table', tableOutput('summary_table'))
                
                ))))

server <- function (input, output) {
  
  output$bar <- renderPlot({
    
    goi = input$entered_gene
    
    #creating the two 'sub' dataframes from the master dataframe based on the gene that was searched
    searched_gene <- df[df['gene'] == goi, ]
    searched_gene_grouped <- df_grouped[df_grouped['gene'] == goi, ]
    
    
    #reordering the x labels such that they are in chronological order
    searched_gene$age <- factor(searched_gene$age, levels = c('E18', 'P4', 'P14', 'P60','P60 + LPS'))
    searched_gene$sex <- factor(searched_gene$sex, levels = c('Male', 'Female'))
    searched_gene_grouped$age <- factor(searched_gene_grouped$age, levels = c('E18', 'P4', 'P14', 'P60','P60 + LPS'))
    searched_gene_grouped$sex <- factor(searched_gene_grouped$sex, levels = c('Male', 'Female'))
 
    ggplot(searched_gene, aes(age, expression, fill = sex))+
      geom_bar(stat = 'identity', data = searched_gene_grouped, position = position_dodge(.9))+
      geom_jitter(position = position_jitterdodge(jitter.width = 0, dodge.width = 0.9), color = 'black', alpha = 0.7, shape = 21)+
      geom_errorbar(data = searched_gene_grouped, aes(ymin = expression-sem, ymax = expression + sem), width = 0.2, position=position_dodge(.9))+
      labs(title = paste(goi, 'Expression'), align = 'center')+
      xlab('Age/Treatment')+
      ylab('Average Expression')+
      theme_classic()+
      theme(plot.title = element_text(hjust = 0.5, size = 18, color = 'black'))
    
    })
  
  output$violin <- renderPlot({
    
    goi = input$entered_gene
    
    #creating the two 'sub' dataframes from the master dataframe based on the gene that was searched
    searched_gene <- df[df['gene'] == goi, ]
    searched_gene_grouped <- df_grouped[df_grouped['gene'] == goi, ]
    
    
    #reordering the x labels such that they are in chronological order
    searched_gene$age <- factor(searched_gene$age, levels = c('E18', 'P4', 'P14', 'P60','P60 + LPS'))
    searched_gene$sex <- factor(searched_gene$sex, levels = c('Male', 'Female'))
    searched_gene_grouped$age <- factor(searched_gene_grouped$age, levels = c('E18', 'P4', 'P14', 'P60','P60 + LPS'))
    searched_gene_grouped$sex <- factor(searched_gene_grouped$sex, levels = c('Male', 'Female'))
    
    ggplot(searched_gene, aes(age, expression, fill = sex))+
      geom_violin(trim = FALSE)+
      geom_jitter(position = position_jitterdodge(jitter.width = 0, dodge.width = 0.9), color = 'black', alpha = 0.9)+
      labs(title = paste(goi, 'Expression'), align = 'center')+
      xlab('Age/Treatment')+
      ylab('Average Expression')+
      theme_classic()+
      theme(plot.title = element_text(hjust = 0.5, size = 18, color = 'black'))
    
  })
  
  output$summary_table <- renderTable({
  goi = input$entered_gene
  
  #creating the two 'sub' dataframes from the master dataframe based on the gene that was searched
  searched_gene <- df[df['gene'] == goi, ]
  searched_gene_grouped <- df_grouped[df_grouped['gene'] == goi, ]
  
  
  #reordering the x labels such that they are in chronological order
  searched_gene$age <- factor(searched_gene$age, levels = c('E18', 'P4', 'P14', 'P60','P60 + LPS'))
  searched_gene$sex <- factor(searched_gene$sex, levels = c('Male', 'Female'))
  searched_gene_grouped$age <- factor(searched_gene_grouped$age, levels = c('E18', 'P4', 'P14', 'P60','P60 + LPS'))
  searched_gene_grouped$sex <- factor(searched_gene_grouped$sex, levels = c('Male', 'Female'))
  
  searched_gene_grouped})
}

shinyApp(ui = ui, server = server)
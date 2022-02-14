library(shiny)
library(shinyalert)
library(shinyvalidate)
library(rsconnect)
library("DT")
library(tidyr)

ui <- tagList(
  fluidPage(
    #includeScript("./www/text.js"),
    tabsetPanel(
      id="maintab",
      tabPanel(
        "Introduction",
        mainPanel(
          h2("Batch correction - What? Why? How?"),
          p("Hello and welcome! We are pleased to be trusted with your data and present to you a visual summary of your dataset.
          
          Batch effect correction is the procedure of removing variability from your data that is not due to your variable of interest. Batch effects are due to technical differences between your samples, such as the type of instrument or even the technician that ran the sample.
          
          For any concerns about or suggestions to improve this report, let the team know at **GroupHeartBioinformaticsSupport@cshs.org**, and we will consider your feedback. We are happy to support your computational needs.")
        )
      ),
      tabPanel(
        "Settings",
        sidebarPanel(
          sliderInput(inputId = "sample_threshold", 
                      label = "1. Sample Threshold", 
                      value = 0.7, min = 0, max = 1),
          sliderInput(inputId = "expgrp_threshold", 
                      label = "2. Experimental group Threshold", 
                      value = 0.5, min = 0, max = 1),
          sliderInput(inputId = "batch_threshold", 
                      label = "3. Batch Threshold", 
                      value = 0.7, min = 0, max = 1),
          fileInput("norm_file", "4. Choose normalized file", accept = ".txt"),
          fileInput("unnorm_file", "5. Choose unnormalized file", accept = ".txt"),
          fileInput("anno_file", "6. Choose annotation file", accept = ".txt"),
          selectInput("cols_of_int", "7. Choose the columns of interest", choices=c(), multiple=TRUE),
          selectInput("samps_for_corr", "8. Choose the \"experimental group:batch\" combo for correlation", choices=c(), multiple=TRUE),
          textInput("output_prefix", label = "9. Pick a prefix for output files", value = "test"),
          textInput("iRT_prot", label = "10. Enter the name of your iRT protein as listed in your data file", value = "iRT_protein")
        ),
        mainPanel(
          h5(strong("Normalized file preview:")),
          dataTableOutput("contents_norm"),
          
          h5(strong("Unnormalized file preview:")),
          dataTableOutput("contents_unnorm"),
          
          h5(strong("Annotation file preview:")),
          dataTableOutput("contents_anno"), 
        ),
      ),
      tabPanel(
        "Log"
      ),
      tabPanel(
        "Results"
      )
    )
  )
)


server <- function(input, output, session) {
  # To validate errors
  iv <- InputValidator$new()
  
  #norm_file<-reactive({
  #  file <- input$norm_file
  #  if (is.null(file)){
  #    dataread<-data.frame(Description="Please upload normalized data file for a preview")
  #  }else{
  #    norm_file_temp <- read.delim(file$datapath, header = TRUE)
  #    if(c("ProteinName", "PeptideSequence", "FragmentIon") %in% names(norm_file_temp)) {
  #      dataread<-read.delim(file$datapath, header = TRUE)
  #    }else {
  #      dataread<-data.frame(Error="!!!ERROR!!! Please upload normalized data file that contains columns/column-names with ProteinName, PeptideSequence and FragmentIon.")      
  #    }
  #  }
  #})
  #output$contents_norm <- renderDataTable(datatable(norm_file(),
  #                          options = list(pageLength = 10)))
  
  ########### Process norm file
  output$contents_norm <- renderDataTable({
    file <- input$norm_file
    
    if (is.null(file)){
      returnValue(data.frame(Description="Please upload normalized file"))
    } else {
      data<-read.delim(file$datapath, header = TRUE)
      if(c("ProteinName", "PeptideSequence", "FragmentIon") %in% names(data)) {
        return(datatable(data, options = list(pageLength = 10)))
      } else{
        shinyalert("Column Error","Please upload normalized data file that contains columns/column-names with ProteinName, PeptideSequence and FragmentIon.",type="error")
        returnValue(data.frame(Error="Please fix column error and reload the file"))
      }
    }
  })
  
  ########### Process unnorm file
  output$contents_unnorm <- renderDataTable({
    file <- input$unnorm_file
    
    if (is.null(file)){
      returnValue(data.frame(Description="Please upload unnormalized file"))
    } else {
      data<-read.delim(file$datapath, header = TRUE)
      if(c("ProteinName", "PeptideSequence", "FragmentIon") %in% names(data)) {
        return(datatable(data, options = list(pageLength = 10)))
      } else{
        shinyalert("Column Error","Please upload unnormalized data file that contains columns/column-names with ProteinName, PeptideSequence and FragmentIon.",type="error")
        returnValue(data.frame(Error="Please fix column error and reload the file"))
      }
    }
  })
  
  ########### Process anno file
  output$contents_anno <- renderDataTable({
    file <- input$anno_file
    
    if (is.null(file)){
      returnValue(data.frame(Description="Please upload annotation file"))
    } else {
      data<-read.delim(file$datapath, header = TRUE)
      if(c("attribute_ExperimentalGroup", "Level3", "order") %in% names(data)) {
        return(datatable(data, options = list(pageLength = 10)))
      } else{
        shinyalert("Column Error","Please upload annotation file that contains columns/column names with attribute_ExperimentalGroup, Level3, and order.",type="error")
        returnValue(data.frame(Error="Please fix column error and reload the file"))
      }
    }
  })
  
  ########### Process anno file to load cols of int and samples for correlation
  anno_file<-reactive({
    file <- input$anno_file
    if (is.null(file)){
      dataread<-data.frame(Description="Please upload annotation file")
    }else{
      anno_file_temp <- read.delim(file$datapath, header = TRUE)
      if(c("attribute_ExperimentalGroup", "Level3", "order") %in% names(anno_file_temp)) {
        dataread<-read.delim(file$datapath, header = TRUE)
      }else {
        dataread<-data.frame(Error="Please fix column error and reload file")
      }
    }
  })
  
  ########### Process cols_of_int
  observe({
    updateSelectInput(session, "cols_of_int",
                      choices = colnames(anno_file()),
                      selected = colnames(anno_file()[0]))
  })
  rv <- reactiveVal(NULL)
  observeEvent(input$cols_of_int, {
    rv(input$cols_of_int)
  })
  
  ########### Process samps_for_corr
  combos <- reactive(NULL)
  observeEvent(input$cols_of_int, {
    anno_file_temp <- anno_file()
    unique_exp_grp <- unique(anno_file_temp$attribute_ExperimentalGroup)
    cols_of_int_temp <- rv()
    combos <- tidyr::expand_grid(unique_exp_grp, cols_of_int_temp)
    combos <- combos %>%
      unite("combo", unique_exp_grp:cols_of_int_temp, sep= ":", 
            remove = FALSE)
    
    updateSelectInput(session, "samps_for_corr",
                      choices = combos$combo,
                      selected = combos$combo[0])
  })
  
  ########### Process output file prefix
  iv$add_rule(
    "output_prefix",
    sv_regex("^[a-zA-Z0-9]*$", "Only alphanumeric characters allowed")
  )
  # `enable()` the validation rules and 
  iv$enable()
  output$prefix_chosen <- renderText({input$output_prefix})
  
}

shinyApp(ui = ui, server = server)
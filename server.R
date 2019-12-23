#
# ########################## #
#
# Loading libraries
#
# ########################## #
# 

# library(graphics)
#library(pheatmap)
library(randomcoloR)
library(shiny)
library(shinyBS)
library(shinysky)
library(shinyjs)
library(V8)
library(shinythemes)
library(plotly)
library(shinyWidgets)
library(shinydashboard)
library(limma)
library(heatmaply)
library(ggplot2)
library(ggsignif)
library(RColorBrewer)
library(randomForest)
library(GenVisR)
library(doBy)
library(parallel)
library(caret)
library(survival)
library(survminer)


function(input, output, session) {
  includeScript("www/ITEMS/js/CNApp-analytics.js")
  useShinyjs()
  options(shiny.maxRequestSize=900*1024^2)
  
  cin_part_run <<- "NO"
  cn_profile_part_run <<- "NO"
  
  all_read_part <- NULL
  all_cin_part <- NULL
  all_profile_part <- NULL
  
  dendro <<- "row"
  
  fixed_colors <<- randomColor(count=1000, hue=c("random"), luminosity=c("bright"))
  # cl <<- makeCluster( ceiling( detectCores()*0.75 ) )
  
  
  observeEvent(input$clear, {
    session$reload()
  })
  observeEvent(input$empty, {
    session$reload()
  })
  
  ###############App fixed terms########
  
  
  ## App fixed terms ##
  dir_funs <<- "funs"
  # directory for functions in app
  
  dir_scripts <- "scripts"
  # directory for scripts in app
  
  demo_dir <- "demo"
  # directory for functions in app
  
  aux_files_dir <- "aux_files"
  # directory for functions in app
  
  demo_data_file <- "demo_data_random160COAD_cms_SURVIVAL_14.11.2018.txt"
  # directory for functions in app
  
  ## App fixed terms ##
  ###########################
  
  #This function is responsible for loading in the selected file
  filedata <- reactive({
    infile <- input$data_browse #file loaded
    if (is.null(infile)) {
      # User has not uploaded a file yet
      return(NULL)
    } else {return(infile)}
  })
  
  
  # # #The following functions prepare output for 'file specifications'
  output$file_specifics <- renderUI({
    # infile <- input$data_browse #file loaded
    if (!is.null(filedata()) | input$use_demo==TRUE) {
      
      
      div(
        h4("File specifications"),
        fluidRow(
          column(3,
                 selectInput(width="100%","select_dec", "Decimal character", 
                             choices = list( "." = ".", "," = ","), selected = ".")
          ),
          column(3,
                 selectInput(width="100%","select_sep", "Column field", 
                             choices = list("tab" = "\t", "," = ","), selected = "\t")
          ),
          column(6,
                 selectInput("select_hg", "Human genome build", width="100%",
                             choices = list("GRCh38/hg38" = "hg38", "GRCh37/hg19" = "hg19"), selected = "hg38")
          )
        ),
        fileInput("annot_data_browse", "Add annotation/clinical data", width="40%",
                  accept = c(
                    "text/csv",
                    "text/comma-separated-values,text/plain",
                    ".csv",
                    ".bed",
                    ".tsv",
                    ".txt/tab-separeted-values"),
                  buttonLabel = "Browse"),
        
        div(style="text-align:left;padding-top:1em;", actionBttn("button_read_data", "Read data", style="simple", size="sm", color="primary")),
        busyIndicator(text="Reading data", wait=200)
      )
      
    }
  })
  
  
  # Reading file data
  observeEvent(input$button_read_data, {
    
    ###############Inputs from user##########
    ## Input from user ##
    decimal_separator <- input$select_dec
    
    field_separator <- input$select_sep
    
    genome_build <- input$select_hg
    # genome build version (options: hg19 (build 37) or hg38)
    
    
    ###
    
    ## Input from user ##
    
    
    
    
    #
    # ################ #
    #
    # Application core
    #
    # ################ #
    #
    
    
    
    ###############Reading data user###################### 
    
    ## Reading data user ##
    
    infile <- input$data_browse #file loaded
    SELECT_DEC <- input$select_dec
    SELECT_SEP <- input$select_sep
    if (!is.null(infile) & input$use_demo==TRUE){
      output$error_demo_read <- renderText({
        paste("Please, deselect 'Demo data' option.", sep="")
      })
      to_return <- NULL
    }else if (is.null(infile) & input$use_demo==FALSE) {
      
    } else {
      
      output$error_demo_read <- renderText({
        paste("", sep="")
      })
      
      if (input$use_demo==TRUE & is.null(infile)){
        df <- read.table(paste(demo_dir, demo_data_file, sep="/"), sep="\t", header=T, dec=".")
      } else {
        df <- read.table(infile$datapath, sep=SELECT_SEP, header=T, dec=SELECT_DEC)
      }
      
      output$error_reading_annot_data <- renderText(
        ""
      )
      
      condition1 <- dim(df)[2]>=5 # 5 our more columns ('ID', 'chr', 'loc.start', 'loc.end' and 'seg.mean' and extra variables)
      to_return <- NULL
      if (condition1=="TRUE"){
        condition2 <- length(grep("ID", colnames(df)))>0 & 
          length(grep("chr", colnames(df)))>0 & 
          length(grep("loc.start", colnames(df)))>0 &
          length(grep("loc.end", colnames(df)))>0 # correct header: we need 'ID', 'chr', 'loc.start', 'loc.end' 
        condition3 <- class(df[,"seg.mean"])=="numeric" # decimal character correct ('logration' column is numeric)
        condition4 <- dim(df[-c(1:5),])[2] <= 26
        ########## UNCOMMENT FOR SERVER VERSION !!!!!!!
        #condition5 <- length(as.character(unique(df[,"ID"]))) <= 160 
      }
      
      if (condition1=="FALSE"){
        output$error_reading_data <- renderText(
          paste("Error reading data. Please, correct your 'file specifications' (column separator = ", SELECT_SEP," ??).",sep="")
        )
        output$user_parameters <- renderUI({})
      } else if (condition2=="FALSE"){ 
        output$error_reading_data <- renderText(
          "Please, check data header. 'ID', 'chr', 'loc.start', 'loc.end' are mandatory terms and must be written correctly."
        )
        output$user_parameters <- renderUI({})
      } else if (condition3=="FALSE"){  
        output$error_reading_data <- renderText(
          paste("Error reading data. Please, correct your 'file specifications' (decimal character = ", SELECT_DEC," ??).",sep="")
        )
        output$user_parameters <- renderUI({})
      } else if (condition4=="FALSE"){  
        output$error_reading_data <- renderText(
          paste("Too much variable columns. Please, reduce them to 26 variable columns top",sep="")
        )
        output$user_parameters <- renderUI({})
      }  
      ########## UNCOMMENT FOR SERVER VERSION !!!!!!!
      # else if (condition2=="TRUE" & condition5=="FALSE"){  
      #   output$error_reading_data <- renderText(
      #    ""
      #   )
      #   output$error_too_many_samples <- renderUI(
      #     HTML("<div style='color:red'>Too much samples for online version of CNApp.<br> Please, reduce number of samples to 160 or download local version at <a href='https://github.com/ait5/CNApp' target='_blank'><b><i>https://github.com/ait5/CNApp</i></b></a>")
      #   )
      #   output$user_parameters <- renderUI({})
      # }  
      # else if (condition1=="TRUE" & condition2=="TRUE" & condition3=="TRUE" & condition4=="TRUE" & condition5=="TRUE"){
      else if (condition1=="TRUE" & condition2=="TRUE" & condition3=="TRUE" & condition4=="TRUE" ){
        output$error_reading_data <- renderText(
          ""
        )
        
        ##Annotation data:
        annot_data_infile <- input$annot_data_browse #file loaded
        SELECT_DEC <- input$select_dec
        SELECT_SEP <- input$select_sep
        if (!is.null(annot_data_infile)){
          annot_df <- read.table(annot_data_infile$datapath, sep=SELECT_SEP, header=T, dec=SELECT_DEC)
          
          condition1 <- dim(annot_df)[2]>=2 # 2 our more columns ('ID', and extra)
          if (condition1=="TRUE"){
            condition2 <- length(grep("ID", colnames(annot_df)))>0 # correct header: we need 'ID' 
            #condition3 <- class(df[,"seg.mean"])=="numeric" # decimal character correct ('logration' column is numeric)
            #condition4 <- dim(df[-c(1:5),])[2] <= 26
          }
          
          if (condition1=="FALSE"){
            output$error_reading_annot_data <- renderText(
              paste("Error reading annotation data. Please, correct your 'file specifications' (column separator = ", SELECT_SEP," ??).",sep="")
            )
          } else if (condition2=="FALSE"){ 
            output$error_reading_annot_data <- renderText(
              "Please, check annotation data header. 'ID', is a mandatory term and must be written correctly."
            )
          } else if (condition1=="TRUE" & condition2=="TRUE"){
            output$error_reading_annot_data <- renderText(
              ""
            )
            
            
            names_in_annot <- as.character(annot_df[,"ID"])
            withProgress({
              n_mat <- as.data.frame(matrix(NA, ncol=ncol(annot_df), nrow=nrow(df)))
              colnames(n_mat) <- colnames(annot_df)
              new_df <- cbind(df, n_mat)
              for (jj in 1:length(names_in_annot)){
                n_sample <- as.character(names_in_annot[jj])
                rows_in_df <- which(df[,"ID"]==n_sample)
                row <- annot_df[jj,]
                for ( xx in 1:ncol(annot_df) ){
                  class_var <- class(annot_df[,xx])
                  term <- as.character(row[,xx])
                  print(term)
                  if (class_var == "numeric" | class_var == "integer"){
                    new_df[rows_in_df, (xx+ncol(df))] <- as.numeric(term)
                  } else if (class_var == "factor" | class_var == "character") {
                    new_df[rows_in_df, (xx+ncol(df))] <- rep(term, length(rows_in_df))
                  }
                  
                }
              }
              df <- new_df[,-(which(colnames(new_df)=="ID")[2])]
            },
            min = 1,
            max= length(names_in_annot),
            value = quantile(1:length(names_in_annot), 0.9),
            message="Adding your annotation data..."
            )
            
          }
        }
        
        output$go_to_funs <- renderUI({
          HTML(
            "<i style=color:gray>Go to app functionalities!&nbsp</i><i class='fa fa-cogs' style='color:orange'></i>"
          )
        })
        ########Data summary#########
        output$summary_data_MAIN <- renderUI(
          HTML("<h4>Summary data:</h4>")
        )
        
        dim_df <- dim(df)
        output$dim_data <- renderUI(
          wellPanel(
            HTML("<h5>Lines: ", dim_df[1], "<br>Columns: ", dim_df[2],"<br><br>N samples: ", length(as.character(unique(df[,"ID"]))), "</h5>")
          )
        )
        
        output$summary_data <- renderPrint(
          summary(df)
        )
        
        
        
        
        ###############Variables from user######################
        
        
        ## Variables from user ##
        data_user <- df
        sample_names <- as.character(unique(data_user[,"ID"])) #obtaining ALL samples names in vector
        
        each_sample <- as.vector(do.call(cbind, mclapply(sample_names, function(x,y){which(y==x)[1]}, y=as.character(data_user[,"ID"]))))
        
        variables_names <- colnames(data_user)[-c(2:5)]
        
        variables_class <- rep(NA, length(variables_names))
        
        list.cont.variables <- list()
        vars_to_eliminate <- c()
        for (g in 1:length(variables_names)){
          var_name <- variables_names[g]
          raw_var <- data_user[each_sample,var_name]
          raw_var[which(raw_var=="")] <- NA
          class_var <- class(raw_var)
          if (class_var=="factor" | class_var=="character"){
            variables_class[g] <- "categoric"
            raw_var_2 <- as.character(raw_var)
            list.cont.variables[[paste(var_name)]] <- raw_var_2	
          } else if (class_var=="numeric" | class_var=="integer") {
            variables_class[g] <- "numeric"
            raw_var_2 <- as.numeric(as.character(raw_var))
            list.cont.variables[[paste(var_name)]] <- raw_var_2
          } else {
            print(paste("Variable ", g, "-", var_name," class can't be coerced...!"))
            vars_to_eliminate <- c(vars_to_eliminate, g)
          }      
        }
        if ( length(vars_to_eliminate)>0 ){
          variables_class <- variables_class[-vars_to_eliminate]
          variables_names <- variables_names[-vars_to_eliminate]
        }
        mat_variables.0 <- as.data.frame(do.call(cbind, list.cont.variables))
        
        if ( ncol(mat_variables.0) == 1 ) {
          mat_variables.0 <- cbind(mat_variables.0, mat_variables.0)
          colnames(mat_variables.0) <- c(colnames(mat_variables.0)[1], "id")
          
          variables_names <- c(variables_names, "id")
          variables_class <- c(variables_class, "categoric")
        }
        
        
        variables_info_mat <- cbind(variables_names, variables_class)
        colnames(variables_info_mat) <- c("name_var", "class_var")
        
        # Categorical variables:
        items_vars_categoric <- which(variables_info_mat[,"class_var"]=="categoric")
        annot_vars_categoric <- variables_info_mat[items_vars_categoric,"name_var"]
        
        for (i in items_vars_categoric){#categoric variables
          mat_variables.0[,i] <- as.character(mat_variables.0[,i])
        }
        
        # Numerical variables:
        items_vars_numeric <- which(variables_info_mat[,"class_var"]=="numeric")
        annot_vars_numeric <- variables_info_mat[items_vars_numeric,"name_var"]
        
        
        for (i in items_vars_numeric){
          mat_variables.0[,i] <- as.numeric(as.character(mat_variables.0[,i]))
        }
        
        
        
        sequential_palettes <- rownames(brewer.pal.info[which(brewer.pal.info[,"category"]=="seq"),])
        diverging_palettes <- rownames(brewer.pal.info[which(brewer.pal.info[,"category"]=="div"),])
        palettes <- c(diverging_palettes, diverging_palettes, diverging_palettes)
        
        variables_info_mat <- cbind(variables_info_mat, rep(NA, nrow(variables_info_mat)))
        colnames(variables_info_mat)[3] <- "color_palette"
        variables_info_mat[,"color_palette"] <- palettes[1:nrow(variables_info_mat)]
        #variables_info_mat[which(variables_info_mat[,"class_var"]=="numeric"),"color_palette"] <- diverging_palettes[1:length(items_vars_numeric)]
        #variables_info_mat[which(variables_info_mat[,"class_var"]=="categoric"),"color_palette"] <- diverging_palettes[1:length(items_vars_categoric)]
        
        ## Variables from user ##
        to_return <- list(df=df, mat_variables.0=mat_variables.0, variables_info_mat=variables_info_mat)
        GLOBAL_DF <<- mat_variables.0
        
        ######Side menus######
        output$side_menu_cin <- renderMenu(
          sidebarMenu(
            menuItem(text = div(HTML("<i class='fa fa-star' style='color:white'></i>&nbsp&nbsp RE-SEG & SCORE &nbsp&nbsp&nbsp&nbsp&nbsp&nbsp<i class='fa fa-cogs' style='color:orange'></i>")), tabName="user_param_cin", selected = TRUE)
          )
        )
        
        output$side_menu_profiling <- renderMenu(
          sidebarMenu(
            menuItem(text = div(HTML("<i class='fa fa-map' style='color:white'></i>&nbsp&nbsp REGION PROFILE &nbsp&nbsp&nbsp&nbsp&nbsp&nbsp<i class='fa fa-cogs' style='color:orange'></i>")), tabName="user_param_seg")
          )
        )
        
        output$side_menu_model <- renderMenu(
          sidebarMenu(
            menuItem(text = div(HTML("<i class='fa fa-microchip' style='color:white'></i>&nbsp&nbsp CLASSIFIER MODEL")), #tabName="home_model"
                     menuSubItem(HTML("<b style=color:orange>CREATE NEW MODEL</b>"), icon=icon("code"), tabName="tab_create_param", selected = F)
                     # menuSubItem(HTML("<b style=color:#00cc66>USE MY MODEL</b>"), icon=icon("user-md"), tabName="tab_my_model_param", selected = F)
            )
          )
        )
        
        ######User parameters #####
        output$user_parameters_cin <- renderUI({
          
          
          # Working with df columns to search for a possible sample_var column
          cols <- colnames(df)
          var_cols <- cols[-c(2:5)]
          possible_names_for_sample_var <- c("sample","Sample","samples","Samples", "ID") # possible names for a sample_var colname
          maybe_sample_var <- as.character(cols[which(var_cols%in%possible_names_for_sample_var)]) # searches for a possible sample_var colname
          if (is.null(maybe_sample_var)){maybe_sample_var <- var_cols[1]}
          
          
          # User parameters panel
          div(
            
            fluidRow(
              column(10,
                     h4("Re-segment data and compute CNA scores FCS, BCS and GCS")
              ),
              column(2,
                     div(style="text-align:right;", actionLink("helpparam_cin","Help", icon=icon("question-circle-o"))),
                     bsModal("modal_param_cin","HELP: User parameters RE-SEG & SCORE","helpparam_cin", includeHTML("./aux_files/help_with_param_cin.html"))
              )
            ),
            
            fluidRow(
              column(4,
                     hidden(div(id="sample-var-re-seg",
                                selectInput("sample_var", "Sample variable",  choices = var_cols, selected=maybe_sample_var) )),
                     
                     br(),
                     
                     radioButtons("def_or_advan", "User parameters:", choices=c("Default", "Advanced"), selected="Default", inline=TRUE)
                     
              )
              #   column(4, 
              #          selectInput("order_by", "Order samples by", choices = var_cols)
              #   )
            )
          )
        })
        
        output$user_parameters_skip_reseg <- renderUI({
          checkboxInput("skip.reseg", label=HTML("<b>Skip re-segmentation</b>"), value=FALSE)
        })
        
        #low.gain.value <- round(log2(2.4/2),2)
        low.gain.value <- 0.2 
        normal.gain.value <- round(log2(3/2),2) 
        high.gain.value <- round(log2(4/2),2) 
        
        
        #low.loss.value <- round(log2(1.6/2),2) 
        low.loss.value <- -0.2
        normal.loss.value <- round(log2(1/2),2) 
        high.loss.value <- round(log2(0.6/2),2) 
        
        output$user_parameters_cin_2 <- renderUI({
          reseg_setts <- "Default"
          reseg_setts <- as.character(input$def_or_advan)
          
          if (reseg_setts=="Default"){
            hidden(div(
              fluidRow(
                
                column(5, 
                       div(id="re-seg_params", wellPanel(
                         helpText("Re-segmentation parameters:"),
                         numericInput("min.length", "Minimal segment length (bp)", value=100000, min = 0, max = Inf, step = 10000, width = NULL),
                         numericInput("max.dist.segm", "Max distance between segments (bp)", value=1000000, min = 0, max = Inf, step = 100000, width = NULL),
                         hidden(div(id="percent-length-to-join",numericInput("percent.dist", "Max distance between segments (in percentage of total length of segments to join)", value=2, min = 0, max = Inf, step = 0.1, width = NULL) )),
                         
                         numericInput("dev.btw.segs", "Max seg.mean deviation between segments", value=0.16, min = 0, max = 1, step = 0.01, width = NULL),
                         numericInput("dev.tozero", "Min seg.mean deviation from zero", value=0.16, min = 0, max = 1, step = 0.01, width = NULL),
                         numericInput("dev.baf", "Max BAF deviation between segments", value=0.1, min = 0, max = 1, step = 0.01, width = NULL)
                       ))
                ),
                column(5,
                       wellPanel(
                         helpText(HTML("CNAs calling thresholds:")),
                         fluidRow(
                           column(6,
                                  numericInput("high.gain", "High gain", value=high.gain.value, step = 0.01, width = NULL),
                                  numericInput("normal.gain", "Normal gain", value=normal.gain.value,  step = 0.01, width = NULL),
                                  numericInput("low.gain", "Low gain", value=low.gain.value,  step = 0.01, width = NULL)
                                  
                           ),
                           column(6,
                                  numericInput("low.loss", "Low loss", value=low.loss.value,  step = 0.01, width = NULL),
                                  numericInput("normal.loss", "Normal loss", value=normal.loss.value, step = 0.01, width = NULL),
                                  numericInput("high.loss", "High loss", value=high.loss.value, step = 0.01, width = NULL)
                           )
                         )
                       ),
                       wellPanel(
                         helpText(HTML("CNAs classification:<br><i>(affecting CNA Scores computation)</i>")),
                         sliderInput("chrom.percent", "Minimal coverage for chromosomal CNAs", 0, 1, value=0.9, step = 0.025, round = FALSE, sep = ",",  ticks = TRUE, animate = FALSE),
                         sliderInput("arm.percent", "Minimal coverage for arm-level CNAs", 0, 1, value=0.5, step = 0.025, round = FALSE, sep = ",",  ticks = TRUE, animate = FALSE)
                       )
                ),
                column(2,
                       hidden(div(id="hidden_cin_params",wellPanel(
                         #helpText("Your parameters to calculate the CNV's are:"),
                         fluidRow(
                           column(6,
                                  numericInput("min.baf", "Minimum BAF to CN-LOH", value=0.2, step = 0.01, width = NULL)
                           ),
                           
                           column(6,
                                  
                                  #sliderInput("chrom.percent", "Minimum coverage for Chromosome alteration", 0, 1, value=0.9, step = 0.025, round = FALSE, sep = ",",  ticks = TRUE, animate = FALSE),
                                  
                                  radioButtons("acrocentric", "Ignore p arms of chromosomes 13, 14, 15, 21, 22", choices = c("Yes","No"), selected = "No", inline = TRUE, width = NULL, choiceNames = NULL, choiceValues = NULL),
                                  
                                  #sliderInput("arm.percent", "Minimum coverage for Arm alteration", 0, 1, value=0.5, step = 0.025, round = FALSE, sep = ",",  ticks = TRUE, animate = FALSE),
                                  
                                  div( title="Low cutoff", sliderInput("focal.percent.low", "Cutoffs to define intensities of Focal alterations", 0, 1, value=0.05, step = 0.025, round = FALSE, sep = ",",  ticks = TRUE, animate = FALSE)),
                                  div( title="Medium cutoff", sliderInput("focal.percent.medium", "", 0, 1, value=0.15, step = 0.025, round = FALSE, sep = ",",  ticks = TRUE, animate = FALSE)),
                                  div( title="High cutoff", sliderInput("focal.percent.high", "", 0, 1, value=0.3, step = 0.025, round = FALSE, sep = ",",  ticks = TRUE, animate = FALSE))
                           )
                         )
                       )))
                       
                )
              )
            ))
          } else if (reseg_setts=="Advanced") {
            div(
              fluidRow(
                
                column(5, 
                       div(id="re-seg_params", wellPanel(
                         helpText("Re-segmentation parameters:"),
                         numericInput("min.length", "Minimal segment length (bp)", value=100000, min = 0, max = Inf, step = 10000, width = NULL),
                         numericInput("max.dist.segm", "Max distance between segments (bp)", value=1000000, min = 0, max = Inf, step = 100000, width = NULL),
                         hidden(div(id="percent-length-to-join",numericInput("percent.dist", "Max distance between segments (in percentage of total length of segments to join)", value=2, min = 0, max = Inf, step = 0.1, width = NULL) )),
                         
                         numericInput("dev.btw.segs", "Max seg.mean deviation between segments", value=0.16, min = 0, max = 1, step = 0.01, width = NULL),
                         numericInput("dev.tozero", "Max seg.mean deviation from zero", value=0.16, min = 0, max = 1, step = 0.01, width = NULL),
                         numericInput("dev.baf", "Max BAF deviation between segments", value=0.1, min = 0, max = 1, step = 0.01, width = NULL)
                       ))
                ),
                column(5,
                       wellPanel(
                         helpText(HTML("CNAs calling thresholds:")),
                         fluidRow(
                           column(6,
                                  numericInput("high.gain", "High gain", value=high.gain.value, step = 0.01, width = NULL),
                                  numericInput("normal.gain", "Normal gain", value=normal.gain.value,  step = 0.01, width = NULL),
                                  numericInput("low.gain", "Low gain", value=low.gain.value,  step = 0.01, width = NULL)
                                  
                           ),
                           column(6,
                                  numericInput("low.loss", "Low loss", value=low.loss.value,  step = 0.01, width = NULL),
                                  numericInput("normal.loss", "Normal loss", value=normal.loss.value, step = 0.01, width = NULL),
                                  numericInput("high.loss", "High loss", value=high.loss.value, step = 0.01, width = NULL)
                           )
                         )
                       ),
                       wellPanel(
                         helpText(HTML("CNAs classification:<br><i>(affecting CNA Scores computation)</i>")),
                         sliderInput("chrom.percent", "Minimal coverage for chromosomal CNAs", 0, 1, value=0.9, step = 0.025, round = FALSE, sep = ",",  ticks = TRUE, animate = FALSE),
                         sliderInput("arm.percent", "Minimal coverage for arm-level CNAs", 0, 1, value=0.5, step = 0.025, round = FALSE, sep = ",",  ticks = TRUE, animate = FALSE)
                       )
                ),
                column(2,
                       hidden(div(id="hidden_cin_params",wellPanel(
                         #helpText("Your parameters to calculate the CNV's are:"),
                         fluidRow(
                           column(6,
                                  numericInput("min.baf", "Minimum BAF to CN-LOH", value=0.2, step = 0.01, width = NULL)
                           ),
                           
                           column(6,
                                  
                                  #sliderInput("chrom.percent", "Minimum coverage for Chromosome alteration", 0, 1, value=0.9, step = 0.025, round = FALSE, sep = ",",  ticks = TRUE, animate = FALSE),
                                  
                                  radioButtons("acrocentric", "Ignore p arms of chromosomes 13, 14, 15, 21, 22", choices = c("Yes","No"), selected = "No", inline = TRUE, width = NULL, choiceNames = NULL, choiceValues = NULL),
                                  
                                  #sliderInput("arm.percent", "Minimum coverage for Arm alteration", 0, 1, value=0.5, step = 0.025, round = FALSE, sep = ",",  ticks = TRUE, animate = FALSE),
                                  
                                  div( title="Low cutoff", sliderInput("focal.percent.low", "Cutoffs to define intensities of Focal alterations", 0, 1, value=0.05, step = 0.025, round = FALSE, sep = ",",  ticks = TRUE, animate = FALSE)),
                                  div( title="Medium cutoff", sliderInput("focal.percent.medium", "", 0, 1, value=0.15, step = 0.025, round = FALSE, sep = ",",  ticks = TRUE, animate = FALSE)),
                                  div( title="High cutoff", sliderInput("focal.percent.high", "", 0, 1, value=0.3, step = 0.025, round = FALSE, sep = ",",  ticks = TRUE, animate = FALSE))
                           )
                         )
                       )))
                       
                )
              )
            )
          }
        })
        
        output$user_parameters_cin_3 <- renderUI({
          div(
            ##Button run cin
            div(style="text-align:left;padding-top:2em;", actionBttn("button_run_cin", "Now, run analysis!", style="fill", size="lg", color="primary")),
            busyIndicator(text="Running", wait=200)
          )
        })
        
        
        output$user_parameters <- renderUI({
          
          
          # Working with df columns to search for a possible sample_var column
          cols <- colnames(df)
          var_cols <- cols[-c(2:5)]
          possible_names_for_sample_var <- c("sample","Sample","samples","Samples", "ID") # possible names for a sample_var colname
          maybe_sample_var <- as.character(cols[which(var_cols%in%possible_names_for_sample_var)]) # searches for a possible sample_var colname
          if (is.null(maybe_sample_var)){maybe_sample_var <- var_cols[1]}
          
          # To know different segmentation names
          aux_files_dir <- "aux_files"
          genome_build <- as.character(input$select_hg)
          name_seg_hg <- paste("segmented_files_", genome_build, sep="")
          name_file <- paste(aux_files_dir, "/", name_seg_hg, ".txt", sep="")
          df_segmented_files <- read.table(name_file, sep="\t", header=T)
          seg_names <- df_segmented_files[,"name"]
          
          # User parameters panel
          div(
            
            
            fluidRow(
              column(10,
                     h4("Generate genome-wide region profiles by custome genome-windows analysis")
              ),
              column(2,
                     div(style="text-align:right;", actionLink("helpparam","Help", icon=icon("question-circle-o"))),
                     bsModal("modal_param","HELP: User parameters REGION PROFILE","helpparam", includeHTML("./aux_files/help_with_param.html"))
              )
            ),
            
            fluidRow(
              
              column(6, 
                     hidden(div(id="sample-var-genome-region",
                                selectInput("sample_var", "Sample variable",  choices = var_cols, selected=maybe_sample_var) )),
                     
                     # selectInput("order_by", "Order samples by", choices = var_cols),
                     br(),
                     radioButtons("segs_radio", "Genome windows by:",
                                  choices = seg_names,
                                  width="100%",
                                  inline = T)
                     # selectInput("method_cor", "Correlation method",  
                     #             choices = list("Pearson" = "pearson", "Kendall" = "kendall", "Spearman" = "spearman"), selected="pearson")
              ),
              
              hidden(div(id="genome-region-ths",
                         column(6, 
                                numericInput("gain_th", "Gain threshold", value=0.23, step=0.01),
                                numericInput("loss_th", "Loss threshold", value=-0.23, step=0.01)
                                
                         ) ))
            )
            
          )
        })
        output$user_parameters_2 <- renderUI({
          div(
            fluidRow(
              column(6,
                     HTML("<p style=color:gray><i>Run first <a style=background-color:white>&nbsp<i class='fa fa-star' style='color:orange'></i>&nbsp <i style=color:black>RE-SEG & SCORE &nbsp</i></a> part for extensive options</i></p>")
              ),
              column(6
                     
              )
            )
          )
        })
        
        output$user_parameters_3 <- renderUI({
          div(
            fluidRow(
              column(6
                     
              ),
              column(6,
                     br(),
                     div(style="text-align:left;padding-top:2em;", actionBttn("button_run_profiling", "Now, run analysis!", style="fill", size="lg", color="primary")),
                     busyIndicator(text="Running", wait=200)
              )
            )
          )
        })
        
        # output$home_model <- renderUI({
        #   
        #   # Home model panel
        #   div(
        #     
        #     
        #     fluidRow(
        #       column(6
        #              
        #       ),
        #       column(6,
        #              div(style="text-align:right;", actionLink("helphome_model","Help", icon=icon("question-circle-o"))),
        #              bsModal("modal_home_model","HELP  Work your model","helphome_model", includeHTML("./aux_files/help_work_your_model.html"))
        #       )
        #     ),
        #     br(),
        #     br(),
        #     div(style="text-align:center", h3("What do you want to do?")),
        #     br(),
        #     br(),
        #     fluidRow(
        #       column(6,
        #              div(style="text-align:center",actionBttn("button_create_model", "Create new model", style="fill", size="md", color="warning")),
        #              busyIndicator(text="Running", wait=200)
        #       ),
        #       
        #       column(6,
        #              div(style="text-align:center",actionBttn("button_my_model", "Use my model", style="fill", size="md", color="success")),
        #              busyIndicator(text="Running", wait=200)
        #              
        #       )
        #     )
        #   )
        # })
        
        output$run_parts_in_create_model_cin <- renderUI({
          div(
            HTML("<p style=color:gray><i>Run first <a style=background-color:white>&nbsp<i class='fa fa-star' style='color:orange'></i>&nbsp <i style=color:black>RE-SEG & SCORE &nbsp</i></a> part for extensive parameters</i></p>")
          )
        })
        output$run_parts_in_create_model_profile <- renderUI({
          div(
            HTML("<p style=color:gray><i>Run first <a style=background-color:white>&nbsp<i class='fa fa-map' style='color:orange'></i>&nbsp <i style=color:black>REGION PROFILE &nbsp</i></a> part for extensive parameters</i></p>")
          )
        })
        output$run_parts_in_create_model_data_loaded <- renderUI({
          div(
            br(),
            
            HTML("You will only be able to select variables from <a style=background-color:white>&nbsp<i class='fa fa-folder' style='color:orange'></i>&nbsp <i style=color:black>LOAD DATA &nbsp</i></a>"),
            
            br(),
            br(),
            div(style="text-align:left", title="click 'Load' to load new parameters", actionBttn("button_param_create_model", "Load", size="sm", color="primary", style="simple"))
          )
        })
        
        
      }
      
      
      
    }# else
    
    all_read_part <<- to_return
    
  })
  
  
  observeEvent(input$button_run_cin, {
    # output$anal_run_cin <- renderUI(
    #   div(style="text-align:right;", p("Check the results!"))
    # )    
    hide("side_menu_cin")
    
    output$side_menu_cin.1 <- renderMenu(
      sidebarMenu(
        menuItem("RE-SEG & SCORE", startExpanded =T, icon=icon("star"),
                 menuSubItem("Parameters", tabName="user_param_cin", icon=icon("cogs")),
                 menuSubItem("Results", tabName="results_cin", icon=icon("signal"), selected=T)
        )
      )
    )
    
    l3<-read.table("./aux_files/cytobands_level3_pq.csv",header=TRUE,sep="\t")
    l4<-read.table("./aux_files/cytobands_level4_chrom.csv",header=TRUE,sep="\t")
    
    ###############Inputs from user##########
    ## Input from user ##
    decimal_separator <- input$select_dec
    
    field_separator <- input$select_sep
    
    genome_build <- input$select_hg
    # genome build version (options: hg19 (build 37) or hg38)
    sample_variable <- input$sample_var
    # column name which has sample names
    
    #order_var_name <- input$order_by
    # selected-by-user 'variable' to order samples
    
    ###
    
    ## Input from user ##
    
    
    
    
    
    ###############Reading data user###################### 
    
    ## Reading data user ##
    infile <- input$data_browse #file loaded
    if (is.null(infile) & input$use_demo==TRUE){
      data_user <- read.table(paste(demo_dir, demo_data_file, sep="/"), sep="\t", header=T, dec=".")
    } else {
      data_user <- read.table(infile$datapath, sep=field_separator, header=T, dec=decimal_separator)
    }
    data_user <- data_user[order(data_user[,sample_variable]),]
    sample_names <- as.character(unique(data_user[,sample_variable])) #obtaining ALL samples names in vector
    
    each_sample <- as.vector(do.call(cbind, mclapply(sample_names, function(x,y){which(y==x)[1]}, y=as.character(data_user[,sample_variable]))))
    
    col_sample_variable <- which(colnames(data_user)==sample_variable)
    #col_order_var_name <- which(colnames(data_user)==order_var_name)
    
    columns_of_interest <- c("BAF", "purity")
    
    extra_cols <- which(colnames(data_user)%in%columns_of_interest)
    
    full <- data_user[,c(1:5,extra_cols)]
    
    ## Reading data user ##
    
    ##Annotation data:
    annot_data_infile <- input$annot_data_browse #file loaded
    SELECT_DEC <- input$select_dec
    SELECT_SEP <- input$select_sep
    if (!is.null(annot_data_infile)){
      annot_df <- read.table(annot_data_infile$datapath, sep=SELECT_SEP, header=T, dec=SELECT_DEC)
      
      names_in_annot <- as.character(annot_df[,sample_variable])
      n_mat <- as.data.frame(matrix(NA, ncol=ncol(annot_df), nrow=nrow(data_user)))
      colnames(n_mat) <- colnames(annot_df)
      new_df <- cbind(data_user, n_mat)
      for (jj in 1:length(names_in_annot)){
        n_sample <- as.character(names_in_annot[jj])
        rows_in_df <- which(data_user[,sample_variable]==n_sample)
        row <- annot_df[jj,]
        for ( xx in 1:ncol(annot_df) ){
          class_var <- class(annot_df[,xx])
          term <- as.character(row[,xx])
          print(term)
          if (class_var == "numeric" | class_var == "integer"){
            new_df[rows_in_df, (xx+ncol(data_user))] <- as.numeric(term)
          } else if (class_var == "factor" | class_var == "character") {
            new_df[rows_in_df, (xx+ncol(data_user))] <- rep(term, length(rows_in_df))
          }
          
        }
      }
      data_user <- new_df[,-(which(colnames(new_df)==sample_variable)[2])]
    }
    
    
    
    ###############Variables from user######################
    
    
    ## Variables from user ##
    
    variables_names <- colnames(data_user)[-c(2:5)]
    
    if (length(which(variables_names%in%columns_of_interest))>0){
      variables_names <- variables_names[-which(variables_names%in%columns_of_interest)]
    }
    
    variables_class <- rep(NA, length(variables_names))
    
    list.cont.variables <- list()
    vars_to_eliminate <- c()
    for (g in 1:length(variables_names)){
      var_name <- variables_names[g]
      raw_var <- data_user[each_sample,var_name]
      raw_var[which(raw_var=="")] <- NA
      class_var <- class(raw_var)
      if (class_var=="factor" | class_var=="character"){
        variables_class[g] <- "categoric"
        raw_var_2 <- as.character(raw_var)
        list.cont.variables[[paste(var_name)]] <- raw_var_2	
      } else if (class_var=="numeric" | class_var=="integer") {
        variables_class[g] <- "numeric"
        raw_var_2 <- as.numeric(as.character(raw_var))
        list.cont.variables[[paste(var_name)]] <- raw_var_2
      } else {
        print(paste("Variable ", g, "-", var_name," class can't be coerced...!"))
        vars_to_eliminate <- c(vars_to_eliminate, g)
      }      
    }
    if ( length(vars_to_eliminate)>0 ){
      variables_class <- variables_class[-vars_to_eliminate]
      variables_names <- variables_names[-vars_to_eliminate]
    }
    mat_variables.0 <- as.data.frame(do.call(cbind, list.cont.variables))
    
    variables_info_mat <- cbind(variables_names, variables_class)
    colnames(variables_info_mat) <- c("name_var", "class_var")
    
    # Categorical variables:
    items_vars_categoric <- which(variables_info_mat[,"class_var"]=="categoric")
    annot_vars_categoric <- variables_info_mat[items_vars_categoric,"name_var"]
    
    for (i in items_vars_categoric){#categoric variables
      mat_variables.0[,i] <- as.character(mat_variables.0[,i])
    }
    
    # Numerical variables:
    items_vars_numeric <- which(variables_info_mat[,"class_var"]=="numeric")
    annot_vars_numeric <- variables_info_mat[items_vars_numeric,"name_var"]
    
    
    for (i in items_vars_numeric){
      mat_variables.0[,i] <- as.numeric(as.character(mat_variables.0[,i]))
    }
    
    
    
    sequential_palettes <- rownames(brewer.pal.info[which(brewer.pal.info[,"category"]=="seq"),])
    diverging_palettes <- rownames(brewer.pal.info[which(brewer.pal.info[,"category"]=="div"),])
    palettes <- c(diverging_palettes, diverging_palettes, diverging_palettes)
    
    variables_info_mat <- cbind(variables_info_mat, rep(NA, nrow(variables_info_mat)))
    colnames(variables_info_mat)[3] <- "color_palette"
    variables_info_mat[,"color_palette"] <- palettes[1:nrow(variables_info_mat)]
    #variables_info_mat[which(variables_info_mat[,"class_var"]=="numeric"),"color_palette"] <- diverging_palettes[1:length(items_vars_numeric)]
    #variables_info_mat[which(variables_info_mat[,"class_var"]=="categoric"),"color_palette"] <- diverging_palettes[1:length(items_vars_categoric)]
    
    
    ###############Variables from user
    
    
    
    
    
    ############MARIA PART################
    l3$length<-l3$end-l3$start
    l4$length<-l4$end-l4$start
    
    
    centro<-l3[1:24*2,"start"]
    names(centro)<-l3[1:24*2,"chr"]
    
    l4$cum<-NA
    
    l4_sort<-l4[order(as.numeric(as.character(l4$chr)),na.last = TRUE),]
    
    
    for (i in 1:nrow(l4_sort)) {
      l4_sort$length<-as.numeric(l4_sort$length)
      l4_sort$cum[i]<-sum(l4_sort[1:i,"length"])
    }
    
    
    centro<-centro[order(as.numeric(as.character(names(centro))),na.last = TRUE)]
    
    
    
    plot_file <- function (orig, chroms=l4_sort, centromers=centro, limits=NA, plot.new=FALSE, col=col) {
      if (length(grep("chr", orig$chr))>0){
        orig$chr<-as.numeric(gsub("chr","",orig$chr))
      } else {
        orig$chr<-as.numeric(orig$chr)
      }
      orig$length<-orig$loc.end-orig$loc.start
      orig<-orig[order(orig$chr),]
      orig<-orig[which(!is.na(orig$chr)),]
      orig$adj.start<-NA
      chroms$chr<-as.character(chroms$chr)
      chroms$to.sum<-chroms$cum-chroms$length
      
      if (plot.new) {	
        plot(0,10, ylim=c(-2,2), xlim=c(0,chroms[nrow(chroms),"cum"]) ,pch=19,cex=0.1 , ylab="seg.mean", xlab="Chromosomes", main=orig$ID[1], xaxt="n", xaxs="i")
        axis(1,at=chroms$to.sum+centro, labels=c(1:22,"X","Y"))
        #text(chroms$to.sum+centro, 2, chroms$chr)
      }
      
      abline(v=chroms$cum)
      abline(v=c(chroms$to.sum+centro),lty=2)
      abline(h=0,lty=1)
      abline(h=c(limits),lty=3)
      
      for (i in 1:nrow(orig)) {
        orig[i,"adj.start"]<-orig[i,"loc.start"]+chroms[which(chroms$chr==orig[i,"chr"]),"to.sum"]
        lines( c(orig[i,"adj.start"], orig[i,"adj.start"]+orig[i,"length"]), c(orig[i,"seg.mean"], orig[i,"seg.mean"]) ,lwd=2, col=col  )
      }
      
      if (!is.null(orig$classified)) {
        color<-rgb(1, 0, 0, alpha=0.1)
        upds<-which(orig$type=="CN-LOH")
        if (length(upds)>=1) {
          for (i in 1:length(upds)) {
            xs<-c( orig$adj.start[upds[i]], orig$adj.start[upds[i]], orig$adj.start[upds[i]]+orig$length[upds[i]], orig$adj.start[upds[i]]+orig$length[upds[i]])
            ys<-c( 2.5, -2.5, -2.5, 2.5)
            polygon(xs, ys, col=color, border=NA)
          }
          
        }
        
      }
      
    }
    
    
    
    
    
    
    #######################################
    # Main script
    #######################################
    
    
    if (input$acrocentric=="Yes") {
      
      acrocentrics<-c(13,14,15,21,22)      
      for (i in 1:length(acrocentrics)){
        l4[which(l4$chr==acrocentrics[i]),"length"]<-l4[which(l4$chr==acrocentrics[i]),"length"] - l3[which(l3$chr==acrocentrics[i] & l3$label=="p"),"end"]
      }
      
    }
    
    
    orig_list<-list()
    
    parameters.used<-data.frame(
      
      purity=NA,
      dev=NA,
      dev.baf=NA,
      
      dist.max=NA,
      percent.dist=NA,
      
      low.cutoff.up=NA,
      medium.cutoff.up=NA,
      high.cutoff.up=NA,
      
      low.cutoff.dw=NA,
      medium.cutoff.dw=NA,
      high.cutoff.dw=NA,
      
      max.tozero=NA,
      min.tozero=NA,
      
      min.baf=NA,
      #max.baf<-NA
      
      chrom.percent=NA,
      arm.percent=NA,
      focal.percent.low=NA,
      focal.percent.medium=NA,
      focal.percent.high=NA)
    
    min.length<-input$min.length
    
    #full<-read.table("example_laia_standard.csv",header=TRUE,sep="\t")
    # inFile<-input$csvfile
    # full<-read.csv(inFile$datapath,header=TRUE,sep="\t")
    
    validate(
      need(try(all(c("ID", "chr", "loc.start", "loc.end", "seg.mean") %in% colnames(full))==TRUE  ), "Needed columns : \"ID\", \"chr\", \"loc.start\", \"loc.end\", \"seg.mean\" not found. ")
    )
    
    
    
    full$chr<-gsub("chr","",as.character(full$chr))
    
    full$chr[which(full$chr=="X")]<-23
    full$chr[which(full$chr=="Y")]<-24
    full$chr<-as.numeric(full$chr)
    
    full$length<-full$loc.end-full$loc.start
    
    #full<-full[which(full$length>min.length),]
    
    
    full$ID<-as.character(full$ID)
    full<-full[with(full, order(full$ID, full$chr, full$loc.start) ),  ]
    
    individual<-list()
    files<-as.character(unique(full$ID))
    
    for (i in 1:length(files)) {
      individual[[i]]<-full[which(full$ID==files[i]),]
    }
    
    
    fullscores_list<-list()
    orig_list<-list()
    
    chromosomal<-vector()
    arm<-vector()
    focal<-vector()
    filt_list<-list()
    
    
    
    ### for each file
    files_to_eliminate <- c()
    for (jj in 1:length(files)) {
      
      print(paste("File:",jj))
      
      file<-individual[[jj]]
      
      file$length<-file$loc.end-file$loc.start
      file$chr<-as.character(file$chr)
      
      # if BAF is not defined, set to 0.5 (not affecting to any computation)
      if (is.null(file$BAF)) {file$BAF<-0.5}
      
      orig_list[[jj]]<-file[,c("ID","chr","loc.start","loc.end","seg.mean","length","BAF")]
      
      
      ####################################
      ########  make parameters  #########
      ####################################
      
      n.loops<-100
      max.dist.segm<-input$max.dist.segm  
      percent.dist<-input$percent.dist
      
      chrom.percent<-input$chrom.percent
      arm.percent<-input$arm.percent
      focal.percent.low<-input$focal.percent.low
      focal.percent.medium<-input$focal.percent.medium
      focal.percent.high<-input$focal.percent.high
      
      low.gain<-input$low.gain
      normal.gain<-input$normal.gain
      high.gain<-input$high.gain
      
      low.loss<-input$low.loss
      normal.loss<-input$normal.loss
      high.loss<-input$high.loss
      
      dev.btw.segs <- input$dev.btw.segs
      
      
      min.baf<-input$min.baf
      #max.baf<-NA
      dev.baf<-input$dev.baf
      
      
      ######### PURITY correction (seg.mean adjustment [re-centralization] )
      
      r.lim <- 0.4 # min purity accepted
      loss.lim <- log2((1/2)*r.lim) # min loss value accepted
      
      if (!is.null(file[1,"purity"])) {
        r<-file[1,"purity"]
        if (r < r.lim) {
          r <- r.lim
        }
        
        v <- file$seg.mean
        
        new.seg.values <- function(n, r){
          inside.log <- ( 2^n + (r-1) ) /r
          if (inside.log <= 2^(loss.lim))  {
            inside.log<- 2^(loss.lim)
          }
          new_value <- log2 ( inside.log )
        }
        
        file$seg.mean <- sapply(v, new.seg.values, r=r, simplify=T)
        
      }
      
      
      ######### CAPPING seg.mean values (loss.lim) when NO purity 
      
      v <- file$seg.mean
      new.seg.values <- function(n){
        if (is.na(n)){n <- 0}
        
        if (n < loss.lim)  {
          new_value <- loss.lim
        } else {
          new_value <- n
        }
      }
      file$seg.mean <- sapply(v, new.seg.values, simplify=T)
      
      
      
      r.cutoff<-1
      
      # low.cutoff.up<-log2((low.gain * r.cutoff + 2 * (1-r.cutoff)) / 2)
      # medium.cutoff.up<-log2((normal.gain * r.cutoff + 2 * (1-r.cutoff)) / 2)
      # high.cutoff.up<-log2((high.gain * r.cutoff + 2 * (1-r.cutoff)) / 2)
      # 
      # low.cutoff.dw<-log2((low.loss * r.cutoff + 2 * (1-r.cutoff)) / 2)
      # medium.cutoff.dw<-log2((normal.loss * r.cutoff + 2 * (1-r.cutoff)) / 2)
      # high.cutoff.dw<-log2((high.loss * r.cutoff + 2 * (1-r.cutoff)) / 2)
      
      low.cutoff.up<-low.gain
      medium.cutoff.up<-normal.gain
      high.cutoff.up<-high.gain
      
      low.cutoff.dw<-low.loss
      medium.cutoff.dw<-normal.loss
      high.cutoff.dw<-high.loss
      
      max.tozero<-input$dev.tozero
      min.tozero<-input$dev.tozero
      
      dev.to.use<- dev.btw.segs
      
      
      ### save specific parameters for each file
      parameters.used[jj,]<-NA
      
      
      #parameters.used$purity[jj]<-r
      parameters.used$dev[jj]<-dev.to.use
      parameters.used$dev.baf[jj]<-dev.baf
      parameters.used$dist.max[jj]<-max.dist.segm
      parameters.used$percent.dist[jj]<-percent.dist
      
      parameters.used$low.cutoff.up[jj]<-low.cutoff.up
      parameters.used$medium.cutoff.up[jj]<-medium.cutoff.up
      parameters.used$high.cutoff.up[jj]<-high.cutoff.up
      parameters.used$low.cutoff.dw[jj]<-low.cutoff.dw
      parameters.used$medium.cutoff.dw[jj]<-medium.cutoff.dw
      parameters.used$high.cutoff.dw[jj]<-high.cutoff.dw
      parameters.used$max.tozero[jj]<-max.tozero
      parameters.used$min.tozero[jj]<-min.tozero
      parameters.used$min.baf[jj]<-min.baf
      
      parameters.used$chrom.percent[jj]<-chrom.percent
      parameters.used$arm.percent[jj]<-arm.percent
      parameters.used$focal.percent.low[jj]<-focal.percent.low
      parameters.used$focal.percent.medium[jj]<-focal.percent.medium
      parameters.used$focal.percent.high[jj]<-focal.percent.high
      
      
      
      
      
      
      
      ######### SKIP RE-SEGMENTATION
      
      skip.resegmentation <- input$skip.reseg 
      
      if (skip.resegmentation == TRUE) {
        filt <- file
        
      } else if (skip.resegmentation == FALSE) {  #skip.reseg == FALSE
        
        
        new.file<-list()
        
        
        ## filter minimum length of segments
        new.file[[1]]<-file[which(file$length > min.length),]
        
        if (nrow(new.file[[1]])==0) {  #when sample is ELIMINATED!!
          files_to_eliminate <- c(files_to_eliminate, jj)
          next
        }
        
        
        
        ## set segments near to 0 (defined by min.tozero and max.tozero), to 0
        new.file[[1]]$seg.mean[ which(new.file[[1]]$seg.mean > min.tozero & new.file[[1]]$seg.mean < max.tozero ) ] <- 0
        
        
        ### remove sex chromosomes
        #new.file[[1]]<-new.file[[1]][which(new.file[[1]]$chr!="chrX" &  new.file[[1]]$chr!="chrY"),]
        new.file[[1]]<-new.file[[1]][which(new.file[[1]]$chr!=24),]
        
        
        zz<-2
        new.file[[zz]]<-data.frame(ID=NA, chr=NA,	loc.start=NA,	loc.end=NA, seg.mean=NA, length=NA, BAF=NA)
        
        
        
        for (zz in 2:n.loops) {
          
          print(zz)
          new.file[[zz]]<-data.frame(ID=NA, chr=NA,	loc.start=NA,	loc.end=NA, seg.mean=NA, length=NA, BAF=NA)
          
          k<-1
          flag<-0
          
          for (i in 1:(nrow(new.file[[zz-1]])-1)) {
            
            if (flag==1) {flag<-0; next}
            s1<-new.file[[zz-1]][i,c("ID","chr","loc.start","loc.end","seg.mean","length","BAF")]
            new.file[[zz]][k,]<-s1
            if (nrow( new.file[[zz-1]])<2) {break}
            s2<-new.file[[zz-1]][i+1,c("ID","chr","loc.start","loc.end","seg.mean","length","BAF")]
            
            if (s1$chr==s2$chr) {
              if (s2$seg.mean<(s1$seg.mean+dev.to.use) & s2$seg.mean>(s1$seg.mean-dev.to.use) & (s2$loc.start-s1$loc.end)<max.dist.segm  & s2$BAF<(s1$BAF+dev.baf) & s2$BAF>(s1$BAF-dev.baf) & (s2$loc.start-s1$loc.end)<(s2$length+s1$length)*percent.dist) {
                new.file[[zz]][k,]<-s1
                new.file[[zz]][k,"loc.end"]<-s2$loc.end
                new.file[[zz]][k,"seg.mean"]<-(s1$seg.mean*s1$length+s2$seg.mean*s2$length)/(s1$length+s2$length)
                if (!is.na(s1$BAF)) {
                  new.file[[zz]][k,"BAF"]<-(s1$BAF*s1$length+s2$BAF*s2$length)/(s1$length+s2$length)
                }
                k<-k+1
                flag<-1
                if (i==(nrow(new.file[[zz-1]])-2)) {
                  new.file[[zz]][k,]<-new.file[[zz-1]][i+2,c("ID","chr","loc.start","loc.end","seg.mean","length","BAF")]; break
                }
                next
              } else {
                if (i==nrow(new.file[[zz-1]])-1) {
                  new.file[[zz]][k+1,]<-s2
                  break
                } 
              } 
            } else {
              if (i==nrow(new.file[[zz-1]])-1) {
                new.file[[zz]][k+1,]<-s2
                break
              }
            }
            k<-k+1
            
          }
          
          if (zz>2) {
            if (nrow(new.file[[zz]])==nrow(new.file[[zz-1]])) {break}
          }
          
        }
        
        
        filt<-new.file[[length(new.file)]]
        
      } #else (skip.resegmentation)
      
      #names(new.file)<-1:length(new.file)
      #library(WriteXLS)
      #WriteXLS(new.file,paste("history",files[jj],".xls",sep="_"))
      
      
      filt$length<-filt$loc.end-filt$loc.start
      
      filt$score<-NA
      filt$classified<-NA
      filt$type<-NA
      filt$intensity<-NA
      filt$weight<-NA
      filt$comments<-NA
      
      N_CIN_jj<-0
      arm_jj<-0
      focal_jj<-0
      
      
      
      
      for (i in 1:nrow(filt)) {
        
        #### remove "chr"
        chr<-gsub("chr","",as.character(filt[i,"chr"]))
        
        
        ### check if it is a chrom-level SCNA
        if (filt[i,"length"] > chrom.percent*l4[which(l4$chr==chr),"length"]  ) {
          
          w<-0
          
          ## classify gains
          if (filt[i,"seg.mean"] >= low.cutoff.up &  filt[i,"seg.mean"] < medium.cutoff.up) w<-1
          if (filt[i,"seg.mean"] >= medium.cutoff.up &  filt[i,"seg.mean"] < high.cutoff.up) w<-2
          if (filt[i,"seg.mean"] >= high.cutoff.up ) w<-3
          
          ## classify losses
          if (filt[i,"seg.mean"] <= (low.cutoff.dw) &  filt[i,"seg.mean"] > (medium.cutoff.dw)) w<-1
          if (filt[i,"seg.mean"] <= (medium.cutoff.dw) &  filt[i,"seg.mean"] > (high.cutoff.dw)) w<-2
          if (filt[i,"seg.mean"] <= (high.cutoff.dw)  ) w<-3
          
          if (!is.na(filt[i,"BAF"])) {
            ## classify upd
            if (filt[i,"seg.mean"] >= (low.cutoff.dw) & filt[i,"seg.mean"] <= (low.cutoff.up) & (filt[i,"BAF"]>(0.5+min.baf) | filt[i,"BAF"]<(0.5-min.baf)) ) {
              w<- 2
              filt$type[i]<-"CN-LOH"
            }
          }
          
          
          
          if (w!=0) {
            filt$classified[i]<-"chromosomal"
            if (is.na(filt$type[i])) {
              if (filt[i,"seg.mean"]>0) filt$type[i]<-"Gain"
              else filt$type[i]<-"Loss"
            }
            filt$score[i]<-w
            filt$intensity[i]<-w
            N_CIN_jj<-N_CIN_jj+w    
          }
          
          # skip to new cnv
          next
          
        }
        
        
        ### check if it is an arm-level SCNA
        arms<-l3[which(l3$chr==chr),]
        centromer<-arms[1,"end"]
        l.p<-(min(c(centromer,filt[i,"loc.end"]))-filt[i,"loc.start"])/arms[1,"length"]
        l.q<-(filt[i,"loc.end"]-(max(c(centromer,filt[i,"loc.start"]))))/arms[2,"length"]
        
        if (l.p> arm.percent & l.q > arm.percent) {
          print("WARNING: both arms significant!")
          filt$comments[i]<-"Both arms significant!"  # but only one is counted
        }
        
        if (l.p > arm.percent) {
          
          w<-0
          
          ## classify gains
          if (filt[i,"seg.mean"] >= low.cutoff.up &  filt[i,"seg.mean"] < medium.cutoff.up) w<-1
          if (filt[i,"seg.mean"] >= medium.cutoff.up &  filt[i,"seg.mean"] < high.cutoff.up) w<-2
          if (filt[i,"seg.mean"] >= high.cutoff.up ) w<-3
          
          ## classify losses
          if (filt[i,"seg.mean"] <= (low.cutoff.dw) &  filt[i,"seg.mean"] > (medium.cutoff.dw)) w<-1
          if (filt[i,"seg.mean"] <= (medium.cutoff.dw) &  filt[i,"seg.mean"] > (high.cutoff.dw)) w<-2
          if (filt[i,"seg.mean"] <= (high.cutoff.dw)  ) w<-3
          
          if (!is.na(filt[i,"BAF"])) {
            ## classify upd
            if (filt[i,"seg.mean"] >= (low.cutoff.dw) & filt[i,"seg.mean"] <= (low.cutoff.up) & (filt[i,"BAF"]>(0.5+min.baf) | filt[i,"BAF"]<(0.5-min.baf)) ) {
              w<- 2
              filt$type[i]<-"CN-LOH"
            }
          }
          
          if (w!=0) {
            filt$classified[i]<-"arm"
            if (is.na(filt$type[i])) {
              if (filt[i,"seg.mean"]>0) filt$type[i]<-"Gain"
              else filt$type[i]<-"Loss"
            }
            filt$score[i]<-w
            filt$intensity[i]<-w
            arm_jj<-arm_jj+w    
          }
          
          
          # skip to new CNV
          next
          
        }
        
        if (l.q > arm.percent) {
          
          w<-0
          
          ## classify gains
          if (filt[i,"seg.mean"] >= low.cutoff.up &  filt[i,"seg.mean"] < medium.cutoff.up) w<-1
          if (filt[i,"seg.mean"] >= medium.cutoff.up &  filt[i,"seg.mean"] < high.cutoff.up) w<-2
          if (filt[i,"seg.mean"] >= high.cutoff.up ) w<-3
          
          ## classify losses
          if (filt[i,"seg.mean"] <= (low.cutoff.dw) &  filt[i,"seg.mean"] > (medium.cutoff.dw)) w<-1
          if (filt[i,"seg.mean"] <= (medium.cutoff.dw) &  filt[i,"seg.mean"] > (high.cutoff.dw)) w<-2
          if (filt[i,"seg.mean"] <= (high.cutoff.dw)  ) w<-3
          
          ## classify upd
          if (!is.na(filt[i,"BAF"])) {
            if (filt[i,"seg.mean"] >= (low.cutoff.dw) & filt[i,"seg.mean"] <= (low.cutoff.up) & (filt[i,"BAF"]>(0.5+min.baf) | filt[i,"BAF"]<(0.5-min.baf)) ) {
              w<- 2
              filt$type[i]<-"CN-LOH"
            }
          }    
          
          if (w!=0) {
            filt$classified[i]<-"arm"
            if (is.na(filt$type[i])) {
              if (filt[i,"seg.mean"]>0) filt$type[i]<-"Gain"
              else filt$type[i]<-"Loss"
            }
            filt$score[i]<-w
            filt$intensity[i]<-w
            arm_jj<-arm_jj+w    
          }
          
          
          # skip to new CNV
          next
          
        }
        
        
        
        ### check if it is a focal-level SCNA
        if (filt[i,"seg.mean"]>low.cutoff.up | filt[i,"seg.mean"]<low.cutoff.dw) {
          
          ww<-0
          
          
          ### identify percentages of arms
          arms<-l3[which(l3$chr==chr),]
          centromer<-arms[1,"end"]
          l.p<-(min(c(centromer,filt[i,"loc.end"]))-filt[i,"loc.start"])/arms[1,"length"]
          l.q<-(filt[i,"loc.end"]-(max(c(centromer,filt[i,"loc.start"]))))/arms[2,"length"]
          l.tot <- max(l.p,l.q)
          
          if (l.p>0 & l.q >0) {
            print("WARNING: focal SCNA including centromer!")
            filt$comments[i]<-"Focal SCNA including centromer"  # only the highest weight is counted
          }
          
          
          ## weight of focal SCNA according to its length
          if (l.tot <= focal.percent.low) ww<-1
          if (l.tot <= focal.percent.medium & l.tot > focal.percent.low ) ww<-2
          if (l.tot <= focal.percent.high & l.tot > focal.percent.medium ) ww<-3
          if (l.tot > focal.percent.high  ) ww<-4
          
          w<-0
          
          ## classify gains
          if (filt[i,"seg.mean"] >= low.cutoff.up &  filt[i,"seg.mean"] < medium.cutoff.up) w<-1
          if (filt[i,"seg.mean"] >= medium.cutoff.up &  filt[i,"seg.mean"] < high.cutoff.up) w<-2
          if (filt[i,"seg.mean"] >= high.cutoff.up  ) w<-3
          
          ## classify losses
          if (filt[i,"seg.mean"] <= (low.cutoff.dw) &  filt[i,"seg.mean"] > (medium.cutoff.dw)) w<-1
          if (filt[i,"seg.mean"] <= (medium.cutoff.dw) &  filt[i,"seg.mean"] > (high.cutoff.dw)) w<-2
          if (filt[i,"seg.mean"] <= (high.cutoff.dw)  ) w<-3
          
          
          if (w!=0) {
            filt$classified[i]<-"focal"
            if (filt[i,"seg.mean"]>0) filt$type[i]<-"Gain"
            else filt$type[i]<-"Loss"
            filt$score[i]<-w*ww
            filt$intensity[i]<-w
            filt$weight[i]<-ww
            
            focal_jj<-focal_jj+w*ww    
          }
          
          
          # skip to new cnv
          next
          
        }
        
        
        
        ## classify CN-LOH
        if (!is.na(filt[i,"BAF"])) {
          
          
          ### identify percentages of arms
          arms<-l3[which(l3$chr==chr),]
          centromer<-arms[1,"end"]
          l.p<-(min(c(centromer,filt[i,"loc.end"]))-filt[i,"loc.start"])/arms[1,"length"]
          l.q<-(filt[i,"loc.end"]-(max(c(centromer,filt[i,"loc.start"]))))/arms[2,"length"]
          l.tot <- max(l.p,l.q)
          
          if (l.p>0 & l.q >0) {
            print("WARNING: focal SCNA including centromer!")
            filt$comments[i]<-"Focal SCNA including centromer"  # only the highest weight is counted
          }
          
          
          ww<-0
          
          ## weight of focal SCNA according to its length
          if (l.tot <= focal.percent.low ) ww<-1
          if (l.tot <= focal.percent.medium & l.tot > focal.percent.low ) ww<-2
          if (l.tot <= focal.percent.high & l.tot > focal.percent.medium ) ww<-3
          if (l.tot > focal.percent.high  ) ww<-4
          
          w<-0
          
          if (filt[i,"seg.mean"] >= (low.cutoff.dw) & filt[i,"seg.mean"] <= (low.cutoff.up) & (filt[i,"BAF"]>(0.5+min.baf) | filt[i,"BAF"]<(0.5-min.baf)) ) {
            w<- 2
          }
          
          
          if (w!=0) {
            filt$classified[i]<-"focal"
            filt$type[i]<-"CN-LOH"
            filt$score[i]<-w*ww
            filt$intensity[i]<-w
            filt$weight[i]<-ww
            focal_jj<-focal_jj+w*ww    
          }
          
          
          # skip to new cnv
          next
          
          
        }  
        
        
      }
      
      chromosomal[jj]<-N_CIN_jj
      arm[jj]<-arm_jj
      focal[jj]<-focal_jj
      
      filt_list[[jj]]<-filt
      
    }
    
    print("Ended")
    
    
    
    if (length(files_to_eliminate) > 0){
      #message
      output$samples_down_in_re_seg <- renderUI({
        HTML(paste("<p style='color:gray'><b>IMPORTANT: </b>Sample/s <u>", paste(files[files_to_eliminate], collapse="; "), "</u> has/have been eliminated.<p>", sep=""))
      })
      
      #files re-segmented
      files <- files[-files_to_eliminate]
      filt_list <- filt_list[-files_to_eliminate]
      
      #scores
      chromosomal <- chromosomal[-files_to_eliminate]
      arm <- arm[-files_to_eliminate]
      focal <- focal[-files_to_eliminate]
      
      #annotations
      mat_variables.0 <- mat_variables.0[-files_to_eliminate,]
    }
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    names(filt_list)<-files
    
    
    FCS <- focal
    BCS <- chromosomal + arm
    
    if (length(which(FCS==0))==length(FCS)){
      norm_FCS <- FCS
    } else {	
      norm_FCS <-(FCS-mean(FCS))/sd(FCS)
    }
    
    if (length(which(BCS==0))==length(BCS)){
      norm_BCS <- BCS
    } else {	
      norm_BCS <-(BCS-mean(BCS))/sd(BCS)
    }
    
    GCS <- norm_FCS + norm_BCS
    
    fullscores<-data.frame(arm, chromosomal, FCS, BCS, GCS)
    rownames(fullscores)<-files
    
    #Scores results
    scores<-data.frame(FCS, BCS, GCS)
    rownames(scores)<-files
    
    
    
    
    to.return<-list(scores=scores, fullscores=fullscores, filt_list=filt_list, orig_list=orig_list, parameters.used=parameters.used, th_gain=low.gain, th_loss=low.gain)
    
    ###############
    clinical.parametric <- NULL
    clinical.nonparametric<-NULL
    clinical.data<- NULL
    
    if (ncol(mat_variables.0)>1){
      mat_variables <- as.data.frame(mat_variables.0[which(mat_variables.0[,1]%in%rownames(scores)),])
      
      if (length(grep(sample_variable, colnames(mat_variables)))>0){
        mat_variables <- mat_variables[,-which(colnames(mat_variables)==sample_variable)]
      }
    }  else { 
      
      mat_variables <- NULL
      clinical.parametric <- NULL
      clinical.nonparametric <- NULL
      
    }
    
    if (!is.null(mat_variables)) {
      clinical <- as.data.frame(mat_variables)
      colnames(clinical) <- colnames(mat_variables.0)[-1]
      rownames(clinical) <- files
      
      results.pval.parametric<-data.frame(matrix(NA,ncol=ncol(scores),nrow=ncol(clinical)))
      colnames(results.pval.parametric)<-colnames(scores)
      rownames(results.pval.parametric)<-colnames(clinical)
      
      results.pval.nonparametric<-data.frame(matrix(NA,ncol=ncol(scores),nrow=ncol(clinical)))
      colnames(results.pval.nonparametric)<-colnames(scores)
      rownames(results.pval.nonparametric)<-colnames(clinical)
      
      ## test scores
      
      for (i in 1:ncol(clinical)) {
        name <- colnames(clinical)[i]
        groups<-clinical[rownames(scores),i]
        
        class_groups <- variables_info_mat[which(variables_info_mat[,"name_var"]==name),"class_var"]
        #class_groups <- class(groups)
        
        
        ### if the clinical variable is categoric
        if ((class_groups=="categoric" | class_groups=="factor") & length(which(!is.na(groups)))>3) {
          
          if (length(unique(groups))==2) {
            for (j in 1:ncol(scores)) {
              ## t.test
              value <- t.test(scores[which(groups==unique(groups)[1]),j], scores[which(groups==unique(groups)[2]),j])$p.val
              if (!is.null(value)){results.pval.parametric[i,j]<- value}
              else if (is.null(value)) {results.pval.parametric[i,j]<- NA}
              
              ## non-parametric t.test: wilcoxon test (here does not compute exact pvalue)
              value <- wilcox.test(scores[which(groups==unique(groups)[1]),j], scores[which(groups==unique(groups)[2]),j],exact=FALSE)$p.val
              if (!is.null(value)){results.pval.nonparametric[i,j]<- value}
              else if (is.null(value)) {results.pval.nonparametric[i,j]<- NA}
              
            }
          }
          
          if (length(unique(groups))>2) {
            for (j in 1:ncol(scores)) {
              ### p-value from ANOVA
              value <- summary(aov(scores[,j]~groups))[[1]][1,5]
              if (!is.null(value)){results.pval.parametric[i,j]<- value}
              else if (is.null(value)) {results.pval.parametric[i,j]<- NA}
              
              ### non-parametric ANOVA: kruskal.test
              value <- kruskal.test(scores[,j]~as.factor(groups))$p.val
              if (!is.null(value)){results.pval.nonparametric[i,j]<- value}
              else if (is.null(value)) {results.pval.nonparametric[i,j]<- NA}
              
            }
            
          }
          
        }
        
        ### if both variables are numerical
        if ((class_groups=="numeric" | class_groups=="integer") & length(which(!is.na(groups)))>3) { 
          for (j in 2:ncol(scores)) {
            ### correlation: Pearson
            value <- cor.test(scores[,j],groups,method="pearson")$p.val
            if (!is.null(value)){results.pval.parametric[i,j]<- value}
            else if (is.null(value)) {results.pval.parametric[i,j]<- NA}
            
            ### correlation: Spearman
            value <- cor.test(scores[,j],groups,method="spearman")$p.val
            if (!is.null(value)){results.pval.nonparametric[i,j]<- value}
            else if (is.null(value)) {results.pval.nonparametric[i,j]<- NA}
            
          }
          
          
        }
        
      }
      
      # to.return[["clinical.parametric"]]<-results.pval.parametric
      # to.return[["clinical.nonparametric"]]<-results.pval.nonparametric
      # to.return[["clinical.data"]]<-clinical[rownames(scores),]
      
      clinical.parametric <-results.pval.parametric
      clinical.nonparametric<-results.pval.nonparametric
      clinical.data<-clinical[rownames(scores),]
    }
    
    all_cin_part <<- to.return
    
    GLOBAL_DF <<- cbind(mat_variables.0, scores)
    
    
    cin_part_run <<- "YES"
    
    
    output$user_parameters_2 <- renderUI({
      div(
        fluidRow(
          column(6,
                 radioButtons("use_cin_ranges", label=HTML("Use new segments from <a style=background-color:white>&nbsp<i class='fa fa-star' style='color:orange'></i>&nbsp <i style=color:gray>RE-SEG & SCORE &nbsp</i></a> part:"), 
                              choices=c("Yes", "No"), 
                              selected="Yes",
                              inline = TRUE)
          ),
          column(6,
                 radioButtons("scores_as_annotations", label=HTML("Add scores from <a style=background-color:white>&nbsp<i class='fa fa-star' style='color:orange'></i>&nbsp <i style=color:gray>RE-SEG & SCORE &nbsp</i></a> as annotation variables:"), 
                              choices=c("Yes", "No"), 
                              selected="Yes",
                              inline = TRUE)
          )
        )
      )
    })
    
    
    output$user_parameters_3 <- renderUI({
      if (input$use_cin_ranges=="No"){
        div(
          fluidRow(
            column(6
                   
            ),
            column(6,
                   br(),
                   div(style="text-align:left;padding-top:2em;", actionBttn("button_run_profiling", "Now, run analysis!", style="fill", size="lg", color="primary")),
                   busyIndicator(text="Running", wait=200)
            )
          )
        )
      } else if (input$use_cin_ranges=="Yes") {
        div(
          fluidRow(
            column(6,
                   radioButtons("cin_ranges", label=HTML("Choose CNAs you want to work with:"), 
                                choices=c("All", "Broad (chromosomal and arm-level)", "Focal"), 
                                selected="All",
                                inline = FALSE)  
            ),
            column(6,
                   br(),
                   div(style="text-align:left;padding-top:2em;", actionBttn("button_run_profiling", "Now, run analysis!", style="fill", size="lg", color="primary")),
                   busyIndicator(text="Running", wait=200)
            )
          )
        )
      }
    })
    
    output$run_parts_in_create_model_cin <- renderUI({
      div(
        radioButtons("cin_param_for_model", label=HTML("Activate <a style=background-color:white>&nbsp<i class='fa fa-star' style='color:orange'></i>&nbsp <i style=color:gray>RE-SEG & SCORE &nbsp</i></a> parameters:"), 
                     choices=c("Yes", "No"), 
                     selected="Yes",
                     inline = TRUE)
      )
    })
    output$run_parts_in_create_model_data_loaded <- renderUI({
      div(
        br(),
        div(style="text-align:left", title="click 'Load' to load new parameters", actionBttn("button_param_create_model", "Load", size="sm", color="primary", style="simple"))
      )
    })
    
    
    
    #    if (input$prediction.models!="No prediction") {
    #      #load(paste("models/",input$prediction.models,".RData",sep=""))
    #      load(paste("models/",input$prediction.models,".RData",sep=""))
    #      library(randomForest)
    #      model<-model.cms
    #      new<-predict(model,fullscores)
    #      prob<-predict(model,fullscores,"prob")
    #      print(new)
    #      to.return[["predictions"]]<-data.frame(ID=names(new),predicted.group=new,prob)
    #    }
    
    
    
    #   return(to.return)
    #   
    #   
    #   
    # })
    
    
    ###################################################
    ################ show results #####################
    ###################################################
    
    
    output$choose_sample_to_plot <- renderUI({
      div(
        selectizeInput("sample_to_plot", "", choices = as.character(files), selected = as.character(files)[1], options = list(maxOptions=3000))
      )
    })
    
    item <- reactive({
      term <- as.character(input$sample_to_plot)
      item <- which(files==term)
      item
    })
    
    output$bef_before <- renderPlot({
      
      i <- as.numeric(as.character(item()))
      
      plot_file(orig_list[[i]], limits=parameters.used[i,c("low.cutoff.up","medium.cutoff.up","high.cutoff.up","low.cutoff.dw","medium.cutoff.dw","high.cutoff.dw")], plot.new=TRUE, col="black")
      #plot_file(filt_list[[i]], limits=parameters.used[i,c("low.cutoff.up","medium.cutoff.up","high.cutoff.up","low.cutoff.dw","medium.cutoff.dw","high.cutoff.dw")], plot.new=FALSE, col="red")      
      #plot_file(filt_list[[i]], limits=parameters.used[i,c("low.cutoff.up","medium.cutoff.up","high.cutoff.up","low.cutoff.dw","medium.cutoff.dw","high.cutoff.dw")], plot.new=TRUE, col="red")      
      
    })
    
    output$bef_after <- renderPlot({
      
      i <- as.numeric(as.character(item()))
      
      #plot_file(orig_list[[i]], limits=parameters.used[i,c("low.cutoff.up","medium.cutoff.up","high.cutoff.up","low.cutoff.dw","medium.cutoff.dw","high.cutoff.dw")], plot.new=TRUE, col="black")
      #plot_file(filt_list[[i]], limits=parameters.used[i,c("low.cutoff.up","medium.cutoff.up","high.cutoff.up","low.cutoff.dw","medium.cutoff.dw","high.cutoff.dw")], plot.new=FALSE, col="red")      
      plot_file(filt_list[[i]], limits=parameters.used[i,c("low.cutoff.up","medium.cutoff.up","high.cutoff.up","low.cutoff.dw","medium.cutoff.dw","high.cutoff.dw")], plot.new=TRUE, col="forestgreen")      
      
    })
    
    # Preparing data file with Re-Seg data:
    all.in.one<-filt_list[[1]]
    
    for (i in 2:length(names(filt_list))) {
      all.in.one[1:nrow(filt_list[[i]])+nrow(all.in.one),]<-filt_list[[i]]
    }

    df_ReSeg <- all.in.one[,c("ID", "chr", "loc.start", "loc.end", "seg.mean", "length", "BAF", "classified", "type")]
     if ( length(unique(df_ReSeg$BAF))==1 & unique(df_ReSeg$BAF)==0.5){
       df_ReSeg <- df_ReSeg[,c("ID", "chr", "loc.start", "loc.end", "seg.mean", "length", "classified", "type")] # prensenting no BAF
    }
    
    ######################################
    
    output$ann_segments_all <- downloadHandler(filename = function(){paste("Re-segmented_samples_CNApp_", Sys.time(), ".tsv", sep="")}, 
                                               content=function (ff) {
                                                 
                                                 df_ReSeg_by_patient <- split(df_ReSeg, df_ReSeg$ID)
                                                  names(df_ReSeg_by_patient) <- unique(df_ReSeg$ID)

                                                  df_ReSeg2 <- as.data.frame(do.call(rbind, lapply(names(df_ReSeg_by_patient), function(x){
                                                    name <- x
                                                    pp <- df_ReSeg_by_patient[[name]]
                                                    row_in_vars <- mat_variables.0[which(mat_variables.0$ID==name),-which(colnames(mat_variables.0)=='ID')]

                                                    new_mat <- as.data.frame(matrix(NA, nrow=nrow(pp), ncol=length(row_in_vars)))
                                                    colnames(new_mat) <- names(row_in_vars)
                                                    for (i in 1:nrow(new_mat)){
                                                      new_mat[i,] <- row_in_vars
                                                    }
                                                    pp2 <- cbind(pp, new_mat)
                                                    pp2
                                                  })))
                                                 print(head(df_ReSeg2))
                                                 write.table(df_ReSeg2,ff,sep="\t",quote=FALSE, row.names=FALSE, na="")                                    
                                               })
    
    ### Freq CN plot
    if (!is.null(mat_variables)){
      output$cnfreq_options <- renderUI({
        radioGroupButtons("cnfreq_options", "", choices=c("All samples", "by variable"), selected="All samples")
      })
    } else if (is.null(mat_variables)){
      output$cnfreq_options <- renderUI({
        radioGroupButtons("cnfreq_options", "", choices=c("All samples"), selected="All samples")
      })
    }
    
    
    sub_var_opt <- eventReactive(as.character(input$cnfreq_options), {
      cnfreq_opt <- as.character(input$cnfreq_options)
      if (cnfreq_opt == "by variable") {
        selectInput("sub_variable_group", "", choices=colnames(mat_variables.0))
      }
    })
    output$sub_variable <- renderUI({
      sub_var_opt()
    })
    
    sub_sub_var_opt <- eventReactive(c(as.character(input$cnfreq_options),as.character(input$sub_variable_group)), {
      cnfreq_opt <- as.character(input$cnfreq_options)
      var <- as.character(input$sub_variable_group)
      col_var <- mat_variables.0[,which(colnames(mat_variables.0)==var)]
      class_var <- class(col_var)
      groups_in_var <- sort(unique(as.character(col_var)))
      if (((class_var != "numeric") | (class_var != "integer")) & cnfreq_opt == "by variable") {
        selectizeInput("vars_to_plot_freq", "", choices=groups_in_var, selected=groups_in_var[1:2], multiple=TRUE)
      } else if (cnfreq_opt == "All samples") {
        HTML("")
      }
      
    })
    output$sub_sub_variable <- renderUI({
      sub_sub_var_opt()
    })
    
    
    
    output$cnfreq_plot_button <- renderUI({
      div(
        ##Button run cin
        div(style="text-align:center;padding-top:2em;", actionBttn("button_cnfreq", "Go!", style="simple", size="sm", color="primary")),
        busyIndicator(text="", wait=200)
      )
    })
    
    observeEvent(input$button_cnfreq, {
      
      f_samples <- eventReactive(c(as.character(input$cnfreq_options),as.character(input$sub_variable_group), as.character(input$vars_to_plot_freq)), {
        cnfreq_opt <- as.character(input$cnfreq_options)
        
        if (cnfreq_opt == "by variable") {
          variable <- as.character(input$sub_variable_group)
          terms_in_var <- input$vars_to_plot_freq
          f_samples <- as.character(mat_variables.0[which(mat_variables.0[,variable]%in%terms_in_var),"ID"])
        } else if (cnfreq_opt == "All samples") {
          f_samples <- as.character(mat_variables.0[,"ID"])
        }
        f_samples
      })
      
      w_samples <- which(names(filt_list)%in%as.character(f_samples()))
      list_segs <- filt_list[w_samples]
      
      # output$prova22 <- renderText({
      #   print(f_samples())
      #   length(list_segs)
      # })
      if (length(list_segs)<=1){
        output$few_samples <- renderUI(HTML("Not enougth samples! Change your variables..."))
      } else if (length(list_segs)<150) {
        output$few_samples <- renderUI(HTML(""))
        withProgress({
          fun_name <- "disjoin.to.FQplot" # function name
          my_fun <- paste(fun_name, ".R", sep="")
          source_fun <- paste(dir_funs, "/", my_fun, sep="")
          source(source_fun)
          # sourcing fun into R
          
          fq_df2 <- disjoin.to.FQplot(list_segs, gain=low.cutoff.up, loss=low.cutoff.dw)
          rownames(fq_df2) <- paste(fq_df2$chr, ":", fq_df2$start, "-", fq_df2$end, sep="")},
          
          min = 1,
          max= length(filt_list),
          value = quantile(1:length(filt_list), 0.9),
          message="Computing recurrent segments...!"
        )
        
      } else { # Region profile analysis by 5Mb
        output$few_samples <- renderUI(HTML(""))
        ## Selecting segmentation file ##
        
        fun_name <- "segmentation.selection" # function name
        my_fun <- paste(fun_name, ".R", sep="") # file function name
        source_fun <- paste(dir_funs, "/", my_fun, sep="")
        source(source_fun) # sourcing fun into R
        
        seg_name_for_fq <- "5Mb"
        
        segmented_file_fq <- segmentation.selection(segmentation_name=seg_name_for_fq, genome_build)
        dimensions_seg_file <- as.character(paste(dim(segmented_file_fq)))
        seg_regions_names <- as.character(segmented_file_fq[,4])
        # apply fun to extract segmentation file('segmented_file')
        
        ## Selecting segmentation file ##
        df <- do.call(rbind, list_segs)
        df <- df[,c("chr", "loc.start", "loc.end", "seg.mean", "ID")]
        
        sp <- split(df, df$ID)
        
        fun_name <- "scoring.segments" # function name
        my_fun <- paste(fun_name, ".R", sep="")
        source_fun <- paste(dir_funs, "/", my_fun, sep="")
        source(source_fun)
        # sourcing fun into R
        
        
        withProgress(
          ANS_result_for_fq <- mclapply(sp, scoring.segments, segmented_file=segmented_file_fq),
          min = 1,
          max= length(sample_names),
          value = quantile(1:length(sample_names), 0.9),
          message=HTML("Too many samples... Perfoming 5Mb genome segmentation! Still working...(aprox. 4 minutes for 350 samples)")
        )
        
        ####### Extracting CNA events:
        fun_name <- "classify.values" # function name
        my_fun <- paste(fun_name, ".R", sep="") # file function name
        source_fun <- paste(dir_funs, "/", my_fun, sep="")
        source(source_fun) # sourcing fun into R
        
        values_prova <- mclapply(ANS_result_for_fq, function(x) as.numeric(as.character(x[,4])))
        
        fun.fun <- function(x){
          as.vector(as.character(lapply(x,classify.values, gain=low.cutoff.up, loss=low.cutoff.dw)))
        }
        data_for_fq <- as.data.frame(do.call(cbind, mclapply(values_prova, fun.fun)))
        
        ## Statitstics CNA counts by region
        fun_name <- "fun.grep" # function name
        my_fun <- paste(fun_name, ".R", sep="") # file function name
        source_fun <- paste(dir_funs, "/", my_fun, sep="")
        source(source_fun) # sourcing fun into R
        
        events_order <- c("normal", "loss", "gain")
        counts_0 <- as.matrix(apply(data_for_fq, 1, fun_grep, z=events_order))
        counts_0_freq <- round((counts_0/length(list_segs) ) *100,2)
        counts <- as.data.frame(cbind(seg_regions_names,t(counts_0_freq)))
        colnames(counts)[1] <- "regions"
        counts_by_regions <- counts
        rownames(counts_by_regions) <- seg_regions_names
        
        fq_df2 <- data.frame(segmented_file_fq[,1:3], counts_by_regions[,c("loss", "normal", "gain")])
      }
      
      
      
      
      ############################################
      
      #Preparing 'p_cna_freq_ReSegs':
      fun_name <- "freq.barplot.GainLoss" # function name
      my_fun <- paste(fun_name, ".R", sep="")
      source_fun <- paste(dir_funs, "/", my_fun, sep="")
      source(source_fun)
      # sourcing fun into R
      
      colors <- c("blue", "beige", "red")
      data <- as.data.frame(fq_df2[,c("loss", "normal", "gain")])
      fq_df3 <- as.data.frame(cbind(rownames(data),fq_df2[,c("loss", "gain")]))
      colnames(fq_df3)[1] <- "Region"
      seg_file <- cbind(fq_df2[,1:3],rownames(data))
      
      withProgress({
        p_cna_freq_ReSegs <- freq.barplot.GainLoss(data, colors, segmented_file=seg_file, bar_gap=0)},
        
        min = 1,
        max= nrow(data),
        value = quantile(1:nrow(data), 0.9),
        message="Plotting your CNAs frequencies...!"
      )
      
      output$freq_cn_filt <- renderPlotly({
        p_cna_freq_ReSegs
      })
      
      output$dw_button_cnfreq <- renderUI(
        div(style="text-align:right",
            downloadButton("download_cnfreq_data", label="TSV"),
            downloadButton("download_cnfreq_plot", label="PNG")
            
        )
      )
      
      output$download_cnfreq_plot <- downloadHandler(
        filename= as.character(paste("CNfreq_plot_reseg_", Sys.time(), ".png", sep="")),
        content=function (file){
          p_cna_freq_ReSegs$x$layout$font$size <- 80
          p_cna_freq_ReSegs$x$layout$xaxis$tickfont$size <- 100
          p_cna_freq_ReSegs$x$layout$margin$t <- 300
          p_cna_freq_ReSegs$x$layout$margin$b <- 700
          p_cna_freq_ReSegs$x$layout$margin$l <- 500
          p_cna_freq_ReSegs$x$layout$margin$r <- 500
          export(p_cna_freq_ReSegs, file = file, vwidth = 7000, vheight = 5000, expand=1)
        }
      )
      output$download_cnfreq_data <- downloadHandler(
        filename= as.character(paste("CNAs_freq_", Sys.time(), ".tsv", sep="")),
        content=function (file){
          write.table(x = fq_df3, file = file, sep = "\t", quote=F, row.names=F)
        }
      )
      
      
    })
    
    
    sub_df_ReSeg <- reactive({
      term <- as.character(input$sample_to_plot)
      df_ReSeg[which(df_ReSeg[,"ID"]==term),]
    })
    
    output$sub_data_ReSeg <- downloadHandler(
      filename = reactive({
        term <- as.character(input$sample_to_plot)
        paste(term, "_post_re-segmentation_", Sys.time(), ".tsv", sep="")
      }),
      content=function (ff) {
        write.table(sub_df_ReSeg(),ff,sep="\t",quote=FALSE, row.names=FALSE, na="")                                    
      })
    
    
    output$bef_after_down_one <- downloadHandler(
      filename = reactive({
        term <- as.character(input$sample_to_plot)
        paste(as.character(term), "_pre&post_re-segmentation_", Sys.time(), ".png", sep="")
      }),
      content = function(file) {
        png(file, width = 8000, height = 7000, res=500)
        par(mfrow=c(2,1))
        i <- as.numeric(as.character(item()))
        
        plot_file(orig_list[[i]], limits=parameters.used[i,c("low.cutoff.up","medium.cutoff.up","high.cutoff.up","low.cutoff.dw","medium.cutoff.dw","high.cutoff.dw")], plot.new=TRUE, col="black")
        #plot_file(filt_list[[i]], limits=parameters.used[i,c("low.cutoff.up","medium.cutoff.up","high.cutoff.up","low.cutoff.dw","medium.cutoff.dw","high.cutoff.dw")], plot.new=FALSE, col="red")      
        plot_file(filt_list[[i]], limits=parameters.used[i,c("low.cutoff.up","medium.cutoff.up","high.cutoff.up","low.cutoff.dw","medium.cutoff.dw","high.cutoff.dw")], plot.new=TRUE, col="forestgreen")      
        
        dev.off()
        
      }
    )
    
    
    output$bef_after_down <- downloadHandler(filename = function() { paste("all_samples_plot_befor-after_resegmentation_", Sys.time(), ".pdf", sep="") }, 
                                             function (ff) {
                                               files<-filt_list
                                               orig<-orig_list
                                               
                                               parameters<-parameters.used
                                               pdf(ff, width=30,height=4)
                                               for (i in 1:length(files)) {
                                                 plot_file(orig[[i]], limits=parameters[i,c("low.cutoff.up","medium.cutoff.up","high.cutoff.up","low.cutoff.dw","medium.cutoff.dw","high.cutoff.dw")], plot.new=TRUE, col="black")
                                                 #plot_file(files[[i]], limits=parameters[i,c("low.cutoff.up","medium.cutoff.up","high.cutoff.up","low.cutoff.dw","medium.cutoff.dw","high.cutoff.dw")], plot.new=FALSE, col="red")      
                                                 plot_file(files[[i]], limits=parameters[i,c("low.cutoff.up","medium.cutoff.up","high.cutoff.up","low.cutoff.dw","medium.cutoff.dw","high.cutoff.dw")], plot.new=TRUE, col="green")      
                                                 
                                               }
                                               dev.off()
                                             })
    
    
    ff.1<-data.frame(ID=rownames(scores),  scores)
    colnames(ff.1) <- c("ID", "FCS", "BCS", "GCS")
    if (!is.null(clinical.data)){
      ff <- cbind(ff.1, clinical.data)
    } else {
      ff <- ff.1
    }
    
    output$scores <- renderDataTable({
      ff
    },
    options = list(scrollX=TRUE, pageLength=5, lengthMenu=c(5,10,15,20), paging=TRUE, searching=TRUE, info=FALSE)
    )
    
    
    output$download_scores <- downloadHandler(
      filename=paste("CNApp_CNA_Scores_&_annotations_", Sys.time(), ".tsv", sep=""),
      content=function (file){
        write.table(x = ff, file = file, sep = "\t", quote=F, row.names=F)
      }
    )
    
    output$dw_buttons_plots_indv_segs <- renderUI({
      div(style="text-align:right;",
          downloadButton("sub_data_ReSeg", label="TSV"),
          downloadButton("bef_after_down_one", label="PNG"),
          downloadButton("bef_after_down", label="PDF (all samples)")
      )
    })
    
    # output$download_scores_raw <- downloadHandler(
    #   filename="Scores_raw.tsv",
    #   content=function (file){
    #     ff<-data.frame(ID=rownames(fullscores),  fullscores)
    #     write.table(x = ff, file = file, sep = "\t", quote=F, row.names=F)
    #   }
    # )
    
    
    
    
    
    ### CNA Scores distribution boxplot:
    if ( ncol(ff)>4  ){
      output$box_scores_options <- renderUI({
        radioGroupButtons("box_scores_options", "", choices=c("Global", "by variable"), selected="Global")
      })
    } else if ( ncol(ff)==4 ){
      output$box_scores_options <- renderUI({
        radioGroupButtons("box_scores_options", "", choices=c("Global"), selected="Global")
      })
    }
    
    
    sub_var_score_opt <- eventReactive(as.character(input$box_scores_options), {
      score_opt <- as.character(input$box_scores_options)
      if (score_opt == "by variable") {
        selectInput("sub_variable_score", "", choices=colnames(ff)[-c(1:4)], selected=colnames(ff)[5])
      }
    })
    output$sub_variable_score <- renderUI({
      sub_var_score_opt()
    })
    
    sub_sub_var_score_opt <- eventReactive(c(as.character(input$box_scores_options),as.character(input$sub_variable_score)), {
      score_opt <- as.character(input$box_scores_options)
      var <- as.character(input$sub_variable_score)
      col_var <- ff[,which(colnames(ff)==var)]
      class_var <- class(col_var)
      groups_in_var <- sort(unique(as.character(col_var)))
      if (((class_var != "numeric") | (class_var != "integer")) & score_opt == "by variable") {
        # if ( score_opt == "by variable") {
        selectizeInput("vars_to_plot_scores", "", choices=groups_in_var, selected=groups_in_var, multiple=TRUE)
      } else if (score_opt == "Global") {
        HTML("")
      }
      
    })
    output$sub_sub_variable_score <- renderUI({
      sub_sub_var_score_opt()
    })
    
    
    
    output$score_plot_button <- renderUI({
      div(
        ##Button run boxplot scores
        div(style="text-align:center;padding-top:2em;", actionBttn("button_score", "Go!", style="simple", size="sm", color="primary")),
        busyIndicator(text="", wait=200)
      )
    })
    
    observeEvent(input$button_score, {
      
      ff_box_score <- eventReactive(c(as.character(input$box_scores_options),as.character(input$sub_variable_score), as.character(input$vars_to_plot_scores)), {
        score_opt <- as.character(input$box_scores_options)
        
        if (score_opt == "by variable") {
          variable <- as.character(input$sub_variable_score)
          columns <- c("BCS", "FCS", "GCS", variable)
          terms_in_var <- input$vars_to_plot_scores
          ff_box_score <- as.data.frame(ff[which(ff[,variable]%in%terms_in_var), columns])
        } else if (score_opt == "Global") {
          columns <- c("BCS", "FCS", "GCS")
          ff_box_score <- as.data.frame(ff[,columns])
        }
        ff_box_score
      })
      
      
      fun_name <- "score.boxplot.by.var" # function name
      my_fun <- paste(fun_name, ".R", sep="") # file function name
      source_fun <- paste(dir_funs, "/", my_fun, sep="")
      source(source_fun) # sourcing fun into R
      
      
      cna_scores <- c("BCS", "FCS", "GCS")
      data <- ff_box_score()
      
      p_scores <- mclapply(cna_scores, score.boxplot.by.var, data=data)
      
      p_bcs <- p_scores[[1]]
      xlab_bcs <- p_bcs$labels$x
      
      p_fcs <- p_scores[[2]]
      xlab_fcs <- p_fcs$labels$x
      
      p_gcs <- p_scores[[3]]
      xlab_gcs <- p_gcs$labels$x
      
      output$box_BCS <- renderPlot({
        p_bcs
      })
      
      output$box_FCS <- renderPlot({
        p_fcs
      })
      
      output$box_GCS <- renderPlot({
        p_gcs
      })
      
      
      
      output$dw_button_BCS_plot <- renderUI(
        div(style="text-align:right",
            downloadButton("download_BCS_plot", label="PNG")
            
        )
      )
      output$dw_button_FCS_plot <- renderUI(
        div(style="text-align:right",
            downloadButton("download_FCS_plot", label="PNG")
            
        )
      )
      output$dw_button_GCS_plot <- renderUI(
        div(style="text-align:right",
            downloadButton("download_GCS_plot", label="PNG")
            
        )
      )
      
      
      output$download_BCS_plot <- downloadHandler(
        filename= as.character(paste("BCS_distribution_plot_by_", as.character(xlab_gcs), "_", Sys.time(), ".png", sep="")),
        content=function (file){
          p_bcs$theme$text$size <- 50
          p_bcs$theme$axis.text.x$size <- 30
          
          png(file, width = 7000, height = 5000, res = 400)
          plot(p_bcs)
          dev.off()
        }
      )
      output$download_FCS_plot <- downloadHandler(
        filename= as.character(paste("FCS_distribution_plot_by_", as.character(xlab_gcs), "_", Sys.time(), ".png", sep="")),
        content=function (file){
          p_fcs$theme$text$size <- 50
          p_fcs$theme$axis.text.x$size <- 30
          
          png(file, width = 7000, height = 5000, res = 400)
          plot(p_fcs)
          dev.off()
        }
      )
      output$download_GCS_plot <- downloadHandler(
        filename= as.character(paste("GCS_distribution_plot_by_", as.character(xlab_gcs), "_", Sys.time(), ".png", sep="")),
        content=function (file){
          p_gcs$theme$text$size <- 50
          p_gcs$theme$axis.text.x$size <- 30
          
          png(file, width = 7000, height = 5000, res = 400)
          plot(p_gcs)
          dev.off()
        }
      )
      
      
      
      
      
    })
    
    
    
    
    
    
    
    
    
    ####Heatmap scores plot####
    output$ht_scores_plot_button <- renderUI({
      div(
        ##Button run cin
        div(style="text-align:center;padding-top:2em;", actionBttn("button_ht_scores", "Plot!", style="simple", size="sm", color="primary")),
        busyIndicator(text="", wait=200)
      )
    })
    
    observeEvent(input$button_ht_scores, {
      # Preparing p_dendro:
      fun_name <- "dendro.Ht.ReSeg" # function name
      my_fun <- paste(fun_name, ".R", sep="") # file function name
      source_fun <- paste(dir_funs, "/", my_fun, sep="")
      source(source_fun) # sourcing fun into R
      
      data <- scores
      colnames(data) <- c("FCS", "BCS", "GCS")
      rownames(data) <- rownames(scores)
      data <- t(data)
      
      mat_vars <- mat_variables.0
      
      dendro <- "row"
      
      ht_scores <- dendro.Ht.ReSeg(data, dendro, mat_vars, variables_info_mat)
      
      output$heatmap_scores <- renderPlotly({
        
        ht_scores
        
      })
      # download buttons:
      output$dw_buttons_ht_scores <- renderUI({
        div(style="text-align:right;",
            downloadButton("PNG_ht_scores", "PNG")
        )
      })
      ####PNG_ht_scores
      output$PNG_ht_scores <- downloadHandler(
        filename = paste("Hetmap_scores_", Sys.time(), ".png", sep=""),
        content = function(file) {
          export(ht_scores, file = file, vwidth = 3000, vheight = 2000, expand=1)
        }
      )
      
    })
    
    if (!is.null(clinical.parametric)){
      
      mm <- clinical.parametric
      
      for (i in 1:nrow(mm)) {
        for (j in 1:ncol(mm)) {
          value <- as.numeric(as.character(mm[i,j]))
          
          if (value <= 0.001 & !is.na(value)) {
            value <- paste(round(value,5), "***", sep=" ")
          } else if (value <= 0.01 & !is.na(value)) {
            value <- paste(round(value,5), "**", sep=" ")
          } else if (value <= 0.05 & !is.na(value)) {
            value <- paste(round(value,5), "*", sep=" ")
          } else if (value > 0.05 & !is.na(value)){
            value <- paste(round(value,5), "(ns)", sep=" ")
          } else {
            value <- value
          }
          
          mm[i,j] <- value
        }
      }
      
      ff_pval_parametric <- data.frame(rownames(clinical.parametric),mm)
      colnames(ff_pval_parametric) <- c("Variable", "FCS", "BCS", "GCS")
      
      output$pval.parametric <- renderDataTable({
        ff_pval_parametric
      },
      options = list(lengthChange=TRUE, paging=TRUE, searching=FALSE, info=FALSE)
      )
      
      output$download_pval.parametric <- downloadHandler(
        filename="Pval-parametric.tsv",
        content=function (file){
          ff_pval_parametric <- data.frame(rownames(clinical.parametric),clinical.parametric)
          colnames(ff_pval_parametric) <- c("Variable", "FCS", "BCS", "GCS")
          write.table(x = ff_pval_parametric, file = file, sep = "\t", quote=F, row.names=F)
        }
      )
      
    }
    
    if (!is.null(clinical.nonparametric)){
      
      mm <- clinical.nonparametric
      
      for (i in 1:nrow(mm)) {
        for (j in 1:ncol(mm)) {
          value <- as.numeric(as.character(mm[i,j]))
          
          if (value <= 0.001 & !is.na(value)) {
            value <- paste(round(value,5), "***", sep=" ")
          } else if (value <= 0.01 & !is.na(value)) {
            value <- paste(round(value,5), "**", sep=" ")
          } else if (value <= 0.05 & !is.na(value)) {
            value <- paste(round(value,5), "*", sep=" ")
          } else if (value > 0.05 & !is.na(value)){
            value <- paste(round(value,5), "(ns)", sep=" ")
          } else {
            value <- value
          }
          
          mm[i,j] <- value
        }
      }
      
      ff_pval_nonparametric <- data.frame(rownames(clinical.parametric),mm)
      colnames(ff_pval_nonparametric) <- c("Variable", "FCS", "BCS", "GCS")
      
      
      output$pval.nonparametric <- renderDataTable({
        ff_pval_nonparametric
      },
      options = list(lengthChange=TRUE, paging=TRUE, searching=FALSE, info=FALSE)
      )
      
      output$download_pval.nonparametric <- downloadHandler(
        filename="Pval-nonparametric.tsv",
        content=function (file){
          
          ff_pval_nonparametric <- data.frame(rownames(clinical.nonparametric),clinical.nonparametric)
          colnames(ff_pval_nonparametric) <- c("Variable", "FCS", "BCS", "GCS")
          write.table(x = ff_pval_nonparametric, file = file, sep = "\t", quote=F, row.names=F)
        }
      )
      
    }
    #################################
    ###### Survival analysis ########
    #################################
    
    var_list <- colnames(ff)
    
    if ( length( grep("surv_status", var_list) ) > 0 & length( grep("surv_time", var_list) ) > 0 ){
      
      var_list_2 <- var_list[-c(which(var_list=="ID"), grep("surv_status", var_list), grep("surv_time", var_list) )]
      
      ## Variable selection (by user)
      output$prepare_surv_analysis <- renderUI({
        div(
          
          fluidRow(
            column(4,
                   #variable for survival status
                   selectizeInput("var_status", "Survival STATUS variable", choices= var_list[grep("surv_status", var_list)], multiple=FALSE, selected= var_list[grep("status", var_list)] ),
                   #variable for survival timing
                   selectizeInput("var_time", "Survival TIME variable", choices= var_list[grep("surv_time", var_list)], multiple=FALSE, selected= var_list[grep("time", var_list)] )
            ),
            column(4,
                   #variable to define sample groups in survival analysis
                   selectizeInput("var_surv_groups", "Variable to DEFINE survival groups", choices= var_list_2, multiple=FALSE )
            ),
            column(4
                   
            )
          )
        )
      })
      
      
      ## Assessing variables
     
      subselect_g_var <- eventReactive(input$var_surv_groups, {
        g_var <- as.character(input$var_surv_groups)
        
        if ( g_var=="FCS"| g_var=="BCS" | g_var=="GCS" ) {
          
          g_values <- as.numeric(as.character(ff[,colnames(ff)==g_var]))
          div(
            sliderInput("subsel_g_var", "Groups to compare", min=min(g_values, na.rm = T), 
                        max=max(g_values, na.rm = T), 
                        value = median(g_values, na.rm = T) ) 
          )
          
        } else {
          class_g_var <- variables_info_mat[which(variables_info_mat[,"name_var"]==g_var),"class_var"]
          
          if (class_g_var=="categoric") {
            l_var <- sort(unique(as.character(mat_variables[,colnames(mat_variables)==g_var])))
            div(
              # HTML(paste("<i>(", length(l_var), " groups)</i>", sep="")),
              # br(),
              checkboxGroupInput("subsel_g_var", "Groups to compare", choices = l_var, selected = l_var)
            )
            
          } else if (class_g_var=="numeric"){
            g_values <- round(as.numeric(as.character(mat_variables[,colnames(mat_variables)==g_var])),3)
            div(
              # HTML(paste("<i>(",length(unique(values)), " groups)</i>", sep="")),
              # br(),
              # sliderInput("slider_group_var", "Groups to compare", min=min(values, na.rm = T), max=max(values, na.rm = T), value = c(min(values, na.rm = T),max(values, na.rm = T))),
              sliderInput("subsel_g_var", "Groups to compare", min=min(g_values, na.rm = T), 
                          max=max(g_values, na.rm = T), 
                          value = median(g_values, na.rm = T) ) 
            )
          }
          
        }
        
        
      })
      
      
      output$sub_g_var <- renderUI({
        div(
          subselect_g_var(),
          
          actionBttn("button_run_surv", "Run!", style="simple", size="sm", color="primary"),
          busyIndicator(text="Running", wait=200)
          
        )
      })
      
      
      ## Variable analysis
      observeEvent(input$button_run_surv, {
        mat_vars <- ff
        
        g_var <- as.character(input$var_surv_groups)
        
        if ( g_var=="FCS"| g_var=="BCS" | g_var=="GCS" ) {
          
          value <- as.numeric(as.character(input$subsel_g_var))
          
          surv_var_status <- as.character(input$var_status)
          surv_var_time <- as.character(input$var_time)
          surv_var_groups <- g_var
          #surv_var_sub_groups <- col_group_var
          
          data <- mat_vars[,c("ID", surv_var_status, surv_var_time, surv_var_groups)]
          colnames(data) <- c("id", "status", "time", "variable")
          data$variable <- as.numeric(as.character(data$variable))
          r_low <- which(data$variable < value)
          r_high <- which(data$variable >= value)
          
          data[r_low,"variable"] <- as.character(paste(g_var, "<", value, sep=""))
          data[r_high,"variable"] <- as.character(paste(g_var, ">=", value, sep=""))
          n_groups <- unique(data$variable)
          
          if (g_var == "BCS") {
            colors <- "deepskyblue"
          } else if (g_var == "FCS") {
            colors <- "darkorange"
          } else if (g_var == "GCS") {  
            colors <- "darkkhaki"
          }
          
        } else {
          
          class_g_var <- variables_info_mat[which(variables_info_mat[,"name_var"]==g_var),"class_var"]
          
          if (class_g_var=="numeric") {
            value <- as.numeric(as.character(input$subsel_g_var))
            # min1 <- min(values1, na.rm = T)
            # max1 <- max(values1, na.rm = T)
            # rows_in_annotation1 <- which(as.numeric(as.character(mat_vars[,g_var]))>=min1 & as.numeric(as.character(mat_vars[,g_var]))<=max1)
            # 
            # 
            # col_group_var <- sort(as.numeric(as.character(mat_vars[rows_in_annotation,g_var])))
            # real_groups <- unique(col_group_var)
          } 
          if (class_g_var=="categoric"){
            real_groups <- input$subsel_g_var
            rows_in_annotation <- which(mat_vars[,g_var]%in%real_groups)
            col_group_var <- as.character(mat_vars[rows_in_annotation,g_var])
          }
          
          surv_var_status <- as.character(input$var_status)
          surv_var_time <- as.character(input$var_time)
          surv_var_groups <- g_var
          surv_var_sub_groups <- col_group_var
          
          data <- mat_vars[,c("ID", surv_var_status, surv_var_time, surv_var_groups)]
          colnames(data) <- c("id", "status", "time", "variable")
          
          data <- data[data$variable %in% surv_var_sub_groups,]
          n_groups <- unique(data$variable)
          
          colors <- as.character(variables_info_mat[which(variables_info_mat[,"name_var"]==surv_var_groups),"color_palette"])
          
        }
        


          if ( length( n_groups ) <= 1 ) {
            
            message <- "<i style=color:red>Not enough groups! (Please, select more groups...)</i>"
            output$not_enough_groups_surv <- renderUI({
              div(
                HTML(message)
              )
            })
            
            output$surv_curves_plot <- renderPlot({
              
            })
            
            output$dw_button_survival_plot <- renderUI(
              div(HTML(""))
            )
            
            output$download_survival_plot <- downloadHandler(
              div(HTML(""))
            )
            
          } else if ( length( n_groups ) > 1 ) {
            
            message <- ""
            output$not_enough_groups_surv <- renderUI({
              div(
                HTML(message)
              )
            })
            
            
            # data <- mat_vars[,c("ID", surv_var_status, surv_var_time, surv_var_groups)]
            # colnames(data) <- c("id", "status", "time", "variable")
            # 
            # data <- data[data$variable %in% surv_var_sub_groups,]
            
            # colors <- c("brown", "cyan")
            
            fun_name <- "survival.plot" # function name
            my_fun <- paste(fun_name, ".R", sep="") # file function name
            source_fun <- paste(dir_funs, "/", my_fun, sep="")
            source(source_fun) # sourcing fun into R
            
            p_survival <- survival.plot(data, colors)
            
            # output$message <- renderUI({
            #   div(
            #     HTML(colors)
            #   )
            # })
            
            output$surv_curves_plot <- renderPlot({
              p_survival
            })
            
            
            output$dw_button_survival_plot <- renderUI(
              div(style="text-align:right",
                  downloadButton("download_survival_plot", label="PNG")
                  
              )
            )
            
            output$download_survival_plot <- downloadHandler(
              filename= as.character(paste("Kaplan-Meier_survival_by_", surv_var_groups, "_", Sys.time(), ".png", sep="")),
              content=function (file){
                
                ggexport(p_survival, filename = file, width=1000, height=900, res=150)
                
              }
            )
            
          } # else if ( length( surv_var_sub_groups ) > 1 )
          
          
        
      })
      
    } #if ( length(grep("surv", mat_variables.0[,1]))>0 )
    else {
      message <- "<i style=color:red>Annotation variables to run survival analysis <u><b>not found</b></u>.<br>
      Survival analysis mandatory column headers are '<u>surv_status</u>' and '<u>surv_time</u>'.<br>
      Please, re-evaluate your data annotation (see auxiliary 'Help' buttons).</i>"
      
      output$prepare_surv_analysis <- renderUI({
        div(
          HTML(message)
        )
      })
    }
    
    
    
    
    #  output$predictions.tab <- renderDataTable({
    #    predictions},
    #    options = list(lengthChange=TRUE, paging=TRUE, searching=TRUE, info=FALSE)
    #  )
    #  
    #  output$download_predictions <- downloadHandler(
    #    filename="Predictions.tsv",
    #    content=function (file){
    #      ff<-predictions
    #      write.table(x = ff, file = file, sep = "\t", quote=F, row.names=F)
    #    }
    #  )
    
    
  })#observeEvent(input$button_run_cin{})
  
  
  
  observeEvent(input$button_run_profiling, {
    
    if (cin_part_run=="YES"){
      q_use_cin_ranges <<- as.character(input$use_cin_ranges)
      q_add_cin_scores <- as.character(input$scores_as_annotations)
      th_gain <- all_cin_part$th_gain
      th_loss <- all_cin_part$th_loss*(-1)
      if (q_use_cin_ranges=="Yes"){
        specific_ranges <<- as.character(input$cin_ranges)
      } else if (q_use_cin_ranges=="No") {
        specific_ranges <<- "All"
      }
    } else if (cin_part_run=="NO") {
      q_use_cin_ranges <<- "No"
      q_add_cin_scores <- "No"
      specific_ranges <<- "All"
      th_gain <- 0.2
      th_loss <- -0.2
    }
    
    
    
    #jj_prova <<- "xdf"
    # output$anal_run <- renderUI(
    #   div(style="text-align:right;", p("Check the results!"))
    # )
    
    hide("side_menu_profiling")
    
    output$side_menu_profiling.1 <- renderMenu(
      
      if (q_use_cin_ranges=="Yes") {
        specific_ranges_side <- specific_ranges
        if (specific_ranges_side=="Broad (chromosomal and arm-level)") { specific_ranges_side <- "Broad" }
        sidebarMenu(
          menuItem("REGION PROFILE", startExpanded =T, icon=icon("map"),
                   menuSubItem(HTML("Parameters &nbsp&nbsp&nbsp<i style=color:green>", specific_ranges_side, "</i>&nbsp&nbsp<i class='fa fa-star' style='color:orange'></i>"), tabName="user_param_seg", icon=icon("cogs")),
                   menuSubItem(HTML("CN profiles &nbsp&nbsp&nbsp&nbsp&nbsp<i style=color:orange>by", segmentation_name, "</i>"), tabName="seg_cna_profile_tab", icon=icon("barcode"), selected=T),
                   menuSubItem("Descriptive regions", tabName="limma_tab", icon=icon("tag"))
          )
        )
      } else if (q_use_cin_ranges=="No") {
        sidebarMenu(
          menuItem("REGION PROFILE", startExpanded =T, icon=icon("map"),
                   menuSubItem("Parameters", tabName="user_param_seg", icon=icon("cogs")),
                   menuSubItem(HTML("CN profiles &nbsp&nbsp&nbsp&nbsp&nbsp<i style=color:orange>by", segmentation_name, "</i>"), tabName="seg_cna_profile_tab", icon=icon("barcode"), selected=T),
                   menuSubItem("Descriptive regions", tabName="limma_tab", icon=icon("tag"))
          )
        )
      }
      
    )
    
    output$cn_profiling_TITLE <- renderUI({
      if (q_use_cin_ranges=="Yes") {
        specific_ranges_side <- specific_ranges
        if (specific_ranges_side=="Broad (chromosomal and arm-level)") { specific_ranges_side <- "Broad" }
        HTML("<h4>  Copy number region profile <i style=color:orange>(by", segmentation_name, ")</i></h4>
           <h5 style=text-align:left;color:black>Using <i><b style=color:green>", specific_ranges_side, " CNAs</b></i> from &nbsp<i class='fa fa-star' style='color:orange'></i>&nbsp <i style=color:gray>RE-SEG & SCORE &nbsp</i></h5>
           <br>")
      } else if (q_use_cin_ranges=="No") {
        HTML("<h4>  Copy number region profile <i style=color:orange>(by", segmentation_name, ")</i></h4>
           <br>
           <h4></h4>")
      }
    })
    
    # output$descriptive_regions_TITLE <- renderUI({
    #   HTML("<h3>Descriptive regions between variable groups <i style=color:orange>(by", segmentation_name, ")</i></h3>")
    # })
    
    output$descriptive_regions_TITLE <- renderUI({
      if (q_use_cin_ranges=="Yes") {
        specific_ranges_side <- specific_ranges
        if (specific_ranges_side=="Broad (chromosomal and arm-level)") { specific_ranges_side <- "Broad" }
        HTML("<h4> Descriptive regions between variable groups <i style=color:orange>(by", segmentation_name, ")</i></h4>
           <h5 style=text-align:left;color:black>Using <i><b style=color:green>", specific_ranges_side, " CNAs</b></i> from &nbsp<i class='fa fa-star' style='color:orange'></i>&nbsp <i style=color:gray>RE-SEG & SCORE &nbsp</i></h5>
           <br>")
      } else if (q_use_cin_ranges=="No") {
        HTML("<h4> Descriptive regions between variable groups <i style=color:orange>(by", segmentation_name, ")</i></h3>
           <br>
           <h4></h4>")
      }
    })
    
    
    
    ###############Inputs from user##########
    ## Input from user ##
    decimal_separator <- input$select_dec
    
    field_separator <- input$select_sep
    
    genome_build <- input$select_hg
    # genome build version (options: hg19 (build 37) or hg38)
    sample_variable <- input$sample_var
    # column name which has sample names
    gain_th <- as.numeric(input$gain_th)
    # gain threshold value selected by user
    loss_th <- as.numeric(input$loss_th)
    # loss threshold value selected by user
    
    # method_cor <- input$method_cor
    # correlation sample method in unsupervised analysis
    # (possibilities: pearson [default, kendall, spearman])
    #order_var_name <- as.character(input$order_by)
    order_var_name <- as.character(input$sample_var)
    # selected-by-user 'variable' to order samples
    segmentation_name <- input$segs_radio
    # segmentation type selected by user to be used as genome segmentation
    
    ###
    
    ## Input from user ##
    
    
    
    
    
    
    
    
    
    
    
    
    #
    # ################ #
    #
    # Application core
    #
    # ################ #
    #
    
    
    
    ###############Reading data user###################### 
    
    ## Reading data user ##
    infile <- input$data_browse #file loaded
    if (is.null(infile) & input$use_demo==TRUE){
      data_user <- read.table(paste(demo_dir, demo_data_file, sep="/"), sep="\t", header=T, dec=".")
    } else {
      data_user <- read.table(infile$datapath, sep=field_separator, header=T, dec=decimal_separator)
    }
    data_user <- data_user[order(data_user[,sample_variable]),]
    sample_names <- as.character(unique(data_user[,sample_variable])) #obtaining ALL samples names in vector
    
    each_sample <- as.vector(do.call(cbind, mclapply(sample_names, function(x,y){which(y==x)[1]}, y=as.character(data_user[,sample_variable]))))
    
    ##Annotation data:
    annot_data_infile <- input$annot_data_browse #file loaded
    SELECT_DEC <- input$select_dec
    SELECT_SEP <- input$select_sep
    if (!is.null(annot_data_infile)){
      annot_df <- read.table(annot_data_infile$datapath, sep=SELECT_SEP, header=T, dec=SELECT_DEC)
      
      names_in_annot <- as.character(annot_df[,sample_variable])
      n_mat <- as.data.frame(matrix(NA, ncol=ncol(annot_df), nrow=nrow(data_user)))
      colnames(n_mat) <- colnames(annot_df)
      new_df <- cbind(data_user, n_mat)
      for (jj in 1:length(names_in_annot)){
        n_sample <- as.character(names_in_annot[jj])
        rows_in_df <- which(data_user[,sample_variable]==n_sample)
        row <- annot_df[jj,]
        for ( xx in 1:ncol(annot_df) ){
          class_var <- class(annot_df[,xx])
          term <- as.character(row[,xx])
          print(term)
          if (class_var == "numeric" | class_var == "integer"){
            new_df[rows_in_df, (xx+ncol(data_user))] <- as.numeric(term)
          } else if (class_var == "factor" | class_var == "character") {
            new_df[rows_in_df, (xx+ncol(data_user))] <- rep(term, length(rows_in_df))
          }
          
        }
      }
      data_user <- new_df[,-(which(colnames(new_df)==sample_variable)[2])]
    }
    
    
    col_sample_variable <- which(colnames(data_user)==sample_variable)
    col_order_var_name <- which(colnames(data_user)==order_var_name)
    
    # re-order columns (ID column to the back! Also order_var_name column...)
    if (col_sample_variable==col_order_var_name){
      data_user <- cbind(data_user[,-col_order_var_name], data_user[,order_var_name])
      colnames(data_user)[length(colnames(data_user))] <- order_var_name
    } else if (col_sample_variable!=col_order_var_name){
      data_user <- cbind(data_user[,-c(col_sample_variable, col_order_var_name)], data_user[,c(col_sample_variable, col_order_var_name)])
    }
    
    
    ## Reading data user ##
    
    
    ###############Variables from user######################
    
    
    ## Variables from user ##
    
    variables_names <- colnames(data_user)[-c(1:4)]
    
    variables_class <- rep(NA, length(variables_names))
    
    list.cont.variables <- list()
    vars_to_eliminate <- c()
    for (g in 1:length(variables_names)){
      var_name <- variables_names[g]
      raw_var <- data_user[each_sample,var_name]
      raw_var[which(raw_var=="")] <- NA
      class_var <- class(raw_var)
      if (class_var=="factor" | class_var=="character"){
        variables_class[g] <- "categoric"
        raw_var_2 <- as.character(raw_var)
        list.cont.variables[[paste(var_name)]] <- raw_var_2	
      } else if (class_var=="numeric" | class_var=="integer") {
        variables_class[g] <- "numeric"
        raw_var_2 <- as.numeric(as.character(raw_var))
        list.cont.variables[[paste(var_name)]] <- raw_var_2
      } else {
        print(paste("Variable ", g, "-", var_name," class can't be coerced...!"))
        vars_to_eliminate <- c(vars_to_eliminate, g)
      }      
    }
    if ( length(vars_to_eliminate)>0 ){
      variables_class <- variables_class[-vars_to_eliminate]
      variables_names <- variables_names[-vars_to_eliminate]
    }
    mat_variables.0 <- as.data.frame(do.call(cbind, list.cont.variables))
    
    if ( ncol(mat_variables.0) == 1 ) {
      mat_variables.0 <- cbind(mat_variables.0, mat_variables.0)
      colnames(mat_variables.0) <- c(colnames(mat_variables.0)[1], "id")
      
      variables_names <- c(variables_names, "id")
      variables_class <- c(variables_class, "categoric")
    }
    
    variables_info_mat <- cbind(variables_names, variables_class)
    colnames(variables_info_mat) <- c("name_var", "class_var")
    
    # Categorical variables:
    items_vars_categoric <- which(variables_info_mat[,"class_var"]=="categoric")
    annot_vars_categoric <- variables_info_mat[items_vars_categoric,"name_var"]
    
    for (i in items_vars_categoric){#categoric variables
      mat_variables.0[,i] <- as.character(mat_variables.0[,i])
    }
    
    # Numerical variables:
    items_vars_numeric <- which(variables_info_mat[,"class_var"]=="numeric")
    annot_vars_numeric <- variables_info_mat[items_vars_numeric,"name_var"]
    
    
    for (i in items_vars_numeric){
      mat_variables.0[,i] <- as.numeric(as.character(mat_variables.0[,i]))
    }
    
    
    
    sequential_palettes <- rownames(brewer.pal.info[which(brewer.pal.info[,"category"]=="seq"),])
    diverging_palettes <- rownames(brewer.pal.info[which(brewer.pal.info[,"category"]=="div"),])
    palettes <- c(diverging_palettes, diverging_palettes, diverging_palettes)
    
    variables_info_mat <- cbind(variables_info_mat, rep(NA, nrow(variables_info_mat)))
    colnames(variables_info_mat)[3] <- "color_palette"
    variables_info_mat[,"color_palette"] <- palettes[1:nrow(variables_info_mat)]
    #variables_info_mat[which(variables_info_mat[,"class_var"]=="numeric"),"color_palette"] <- diverging_palettes[1:length(items_vars_numeric)]
    #variables_info_mat[which(variables_info_mat[,"class_var"]=="categoric"),"color_palette"] <- diverging_palettes[1:length(items_vars_categoric)]
    
    ## Variables from user ##
    
    ###############Selecting segmentation file######################
    
    ## Selecting segmentation file ##
    
    fun_name <- "segmentation.selection" # function name
    my_fun <- paste(fun_name, ".R", sep="") # file function name
    source_fun <- paste(dir_funs, "/", my_fun, sep="")
    source(source_fun) # sourcing fun into R
    
    segmented_file <- segmentation.selection(segmentation_name, genome_build)
    dimensions_seg_file <- as.character(paste(dim(segmented_file)))
    seg_regions_names <- as.character(segmented_file[,4])
    # apply fun to extract segmentation file('segmented_file')
    
    ## Selecting segmentation file ##
    
    ###############Use CIN results??######################
    
    ## Use CIN results?? ##
    
    if (q_use_cin_ranges=="Yes") {
      
      if (specific_ranges!="All"){
        if (specific_ranges=="Broad (chromosomal and arm-level)"){
          class_ranges <- c("chromosomal", "arm")
        } else if (specific_ranges=="Focal") {
          class_ranges <- c("focal")
        }
        
        filt_list_1 <- mclapply(all_cin_part[["filt_list"]], function(x){
          
          items <-  which(x[,"classified"]%in%class_ranges)
          if (length(items)>0){
            new_x <-  x[items,]
            new_x[,1:5]
          } else {
            new_x <- as.data.frame(cbind(x[1,1:4], 0))
            colnames(new_x) <- colnames(x)[1:5]
            new_x
          }
          
        })
        filt_list_2 <- filt_list_1
        names_in_cin <- names(filt_list_2)
        
        
        
        
      } else if (specific_ranges=="All") {
        filt_list_2 <- mclapply(all_cin_part[["filt_list"]],function(x){x[,1:5]})
        names_in_cin <- names(filt_list_2)
        
      }
      
      
      rows_to_select <- which(as.character(mat_variables.0[,sample_variable])%in%names_in_cin)
      mat_variables.0 <- mat_variables.0[rows_to_select,]
      
      if (q_add_cin_scores=="Yes"){
        cin_scores <- round(all_cin_part[["scores"]][names_in_cin,],3)
        mat_variables.0 <- cbind(mat_variables.0, cin_scores)
        
        new_mat <- matrix(ncol=3, nrow=ncol(cin_scores))
        new_mat[,1] <- colnames(cin_scores)
        new_mat[,2] <- rep("numeric", length(ncol(cin_scores)))
        new_mat[,3] <- diverging_palettes[1:ncol(cin_scores)]
        variables_info_mat <- rbind(variables_info_mat, new_mat)
      }
      
      
      data_user_2 <- as.data.frame(do.call(rbind, filt_list_2))
      sample_names <- as.character(unique(data_user_2[,sample_variable]))
      
      column_sample_variable <- which(colnames(data_user_2)==sample_variable)
      cols_selected <- c(2:5, column_sample_variable)
      df <- data_user_2[,cols_selected] # chromosome-start-end-value-sample matrix from user input
      
    } else if (q_use_cin_ranges=="No") {
      
      if (q_add_cin_scores=="Yes"){
        cin_scores <- round(all_cin_part[["scores"]],3)
        names_in_cin <- rownames(cin_scores)
        yes_in_cin <- which(sample_names%in%names_in_cin==TRUE)
        not_in_cin <- which(sample_names%in%names_in_cin==FALSE)
        
        if (length(not_in_cin)>0){
          new_mat <- matrix(nrow=dim(mat_variables.0)[1], ncol=ncol(cin_scores))
          colnames(new_mat) <- colnames(cin_scores)
          
          mat_variables.0 <- cbind(mat_variables.0, new_mat)
          
          e_vector <- rep(NA, ncol(cin_scores))			
          mat_variables.0[not_in_cin,colnames(new_mat)] <- e_vector
          
          mat_variables.0[yes_in_cin,colnames(new_mat)] <- cin_scores
          
        } else if (length(not_in_cin)==0) {
          mat_variables.0 <- cbind(mat_variables.0, cin_scores)
        }
        
        
        
        new_mat <- matrix(ncol=3, nrow=ncol(cin_scores))
        new_mat[,1] <- colnames(cin_scores)
        new_mat[,2] <- rep("numeric", length(ncol(cin_scores)))
        new_mat[,3] <- diverging_palettes[1:ncol(cin_scores)]
        variables_info_mat <- rbind(variables_info_mat, new_mat)
      }
      
      
      column_sample_variable <- which(colnames(data_user)==sample_variable)
      cols_selected <- c(1:4, column_sample_variable)
      df <- data_user[,cols_selected] # chromosome-start-end-value-sample matrix from user input
      
      
    }
    
    
    ## Use CIN results?? ##
    
    ###############Applying segmentation assessment and CNA profile generation######################
    
    ## Applying segmentation assessment and CNA profile generation ##
    sp <- split(df, df$ID)
    
    fun_name <- "scoring.segments" # function name
    my_fun <- paste(fun_name, ".R", sep="")
    source_fun <- paste(dir_funs, "/", my_fun, sep="")
    source(source_fun)
    # sourcing fun into R
    
    column_sample_variable <- which(colnames(data_user)==sample_variable)
    cols_selected <- c(1:4, column_sample_variable)
    
    
    
    withProgress(
      ANS_result <- mclapply(sp, scoring.segments, segmented_file),
      min = 1,
      max= length(sample_names),
      value = quantile(1:length(sample_names), 0.9),
      message="Extensive region analyisis! Still working..."
    )
    
    
    ## Applying segmentation assessment and CNA profile generation ##
    
    ###############Extracting 'mean.region' matrix [regions, samples]######################
    
    ## Extracting 'mean.region' matrix [regions, samples] ##
    
    x <- ANS_result
    col_selected <- 4 # 'mean.reg.sample' column
    
    data.0.means <- as.data.frame(do.call(cbind,mclapply(x, function(x) x[,col_selected])))#creating data frame with values of certain column ('col_selected')
    rownames(data.0.means) <- as.character(segmented_file[,4])
    
    GLOBAL_DF <<- cbind(GLOBAL_DF, t(data.0.means))
    ## Extracting 'mean.region' matrix [regions, samples] ##
    
    
    ###############Ordering by selected-by-user 'variable'###################### 
    
    ## Ordering by selected-by-user 'variable' ##
    
    order_by_variable <- order(mat_variables.0[,order_var_name])
    
    mat_variables <- mat_variables.0[order_by_variable,]
    colnames(mat_variables) <- colnames(mat_variables.0)
    
    data.means <- data.0.means[,order_by_variable]
    rownames(data.means) <- as.character(segmented_file[,4])
    
    sample_names <- mat_variables[,sample_variable]
    
    ## Ordering by selected-by-user variable ##
    
    profile_part <- list( sample_variable=sample_variable, variables_info_mat=variables_info_mat, order_by_variable=order_by_variable, segmented_file=segmented_file, mat_variables=mat_variables, data.means=data.means, sample_names=sample_names)
    
    ############### COPY NUMBER PROFILING part run??###################### 
    
    cn_profile_part_run <<- "YES"
    
    ########SIDEBAR with PARAMETERS (plots parameters)#########
    
    output$user_friendly_plot_param_1 <- renderUI({
      div(
        h4("Plot options:"),
        
        selectInput("re_order_by", "Order data by:", choices=as.character(variables_info_mat[,1]), selected = order_var_name),
        
        numericInput("re_gain", "Gain threshold:", value=th_gain, step=0.01),
        numericInput("re_loss", "Loss threshold:", value=th_loss, step=0.01),
        
        selectInput("re_method_correlation", "Correlation method",  
                    choices = list("Pearson" = "pearson", "Kendall" = "kendall", "Spearman" = "spearman"), selected="Pearson"),
        
        hr()
      )
    })
    
    output$user_friendly_plot_param_2 <- renderUI({
      div(
        selectizeInput("tracks", "Annotation track/s:", choices=as.character(variables_info_mat[,1]), selected = c(as.character(variables_info_mat[1,1]), order_var_name), multiple=T)
      )
    })
    
    
    output$user_friendly_plot_param_3 <- renderUI({
      div(
        actionBttn("button_plots_cn_profiles", "Plot!", style="simple", size="sm", color="primary"),
        busyIndicator(text="", wait=200)
      )
    })
    
    
    ########SIDEBAR in LIMMA #########
    
    
    
    output$side_param_limma <- renderUI({
      div(
        h4("Analysis parameters:"),
        selectInput("group_var", "Define groups by",  choices = colnames(mat_variables)[-which(colnames(mat_variables)==sample_variable)])
      )
    })
    
    subselect_group_var <- eventReactive(input$group_var, {
      group_variable <- input$group_var
      class_var <- as.character(variables_info_mat[which(variables_info_mat[,1]==group_variable),2])
      if (class_var=="categoric") {
        legend_var <- sort(unique(as.character(mat_variables[,colnames(mat_variables)==group_variable])))
        div(
          HTML(paste("<i>(", length(legend_var), " groups)</i>", sep="")),
          br(),
          checkboxGroupInput("check_group_var", "Groups to compare", choices = legend_var, selected = legend_var)
        )
      } else if (class_var=="numeric"){
        values <- round(as.numeric(as.character(mat_variables[,colnames(mat_variables)==group_variable])),3)
        id <- paste("slider_", var_name, sep="")
        div(
          HTML(paste("<i>(",length(unique(values)), " groups)</i>", sep="")),
          br(),
          sliderInput("slider_group_var", "Groups to compare", min=min(values, na.rm = T), max=max(values, na.rm = T), value = c(min(values, na.rm = T),max(values, na.rm = T))),
          sliderInput("slider_group_var_2", "", min=min(values, na.rm = T), max=max(values, na.rm = T), value = c(min(values, na.rm = T),max(values, na.rm = T))) 
        )
      }
      
    })
    
    
    output$side_param_limma_2 <- renderUI({
      div(
        subselect_group_var(),
        hr(),
        selectInput("pval_type", "Select p values", choices = c("Adj.p-value", "p-value"), selected = "p-value"),
        sliderInput("slider_pval_limma", "", value=c("0","0.1"), min=0, max=1, step=0.001, round=TRUE),
        hr(),
        actionBttn("button_run_limma", "Run!", style="simple", size="sm", color="primary"),
        busyIndicator(text="Running", wait=200),
        hr(),
        selectizeInput("reg_selected", "Select region to plot", choices = seg_regions_names, options = list(maxOptions=3000))
      )
    })
    
    
    
    
    
    ######WHEN PLOT CN profiles#######
    observeEvent(input$button_plots_cn_profiles, {
      
      output$prova <- renderText({
        # paste("Parameters", field_separator, decimal_character, segmentation_name, genome_build, 
        #       sample_variable, gain_th, loss_th, method_cor, order_var_name, sep=" // ")
        dim(data_user)
        
      })  
      
      
      re_order_var <- as.character(input$re_order_by)
      
      re_gain_th <- as.numeric(as.character(input$re_gain))
      re_loss_th <- as.numeric(as.character(input$re_loss))
      
      var_tracks <- as.character(input$tracks)
      
      re_method_cor <- as.character(input$re_method_correlation)
      
      if (length(var_tracks)>6){
        output$message_text <- renderUI({
          HTML("<i style=color:red>(Max. tracks=6). Please, reduce annotation tracks.</i>")
        })
        
      } else if (length(var_tracks)<1){
        output$message_text <- renderUI({
          HTML("<i style=color:red>Please, select at least one annotation track.</i>")
        })
      } else {
        output$message_text <- renderUI({
          HTML("")
        })
        ## RE-Ordering by selected-by-user 'variable' ##
        
        order_by_variable <- order(mat_variables.0[,re_order_var])
        
        mat_variables <- mat_variables.0[order_by_variable,]
        
        data.means <- data.0.means[,order_by_variable]
        
        sample_names <- mat_variables[,sample_variable]
        
        
        ## RE-Ordering by selected-by-user variable ##
        
        
        ####### Extracting CNA events:
        fun_name <- "classify.values" # function name
        my_fun <- paste(fun_name, ".R", sep="") # file function name
        source_fun <- paste(dir_funs, "/", my_fun, sep="")
        source(source_fun) # sourcing fun into R
        
        values_prova <-mclapply(ANS_result, function(x) as.numeric(as.character(x[,4])))
        
        fun.fun <- function(x){
          as.vector(as.character(lapply(x,classify.values, gain=re_gain_th, loss=re_loss_th)))
        }
        data.0.events <- as.data.frame(do.call(cbind,mclapply(values_prova, fun.fun)))
        data.events <- data.0.events[,order_by_variable]
        ## CNA events plot summary ##
        
        fun_name <- "fun.grep" # function name
        my_fun <- paste(fun_name, ".R", sep="") # file function name
        source_fun <- paste(dir_funs, "/", my_fun, sep="")
        source(source_fun) # sourcing fun into R
        
        
        fun_name <- "extract.stats" # function name
        my_fun <- paste(fun_name, ".R", sep="") # file function name
        source_fun <- paste(dir_funs, "/", my_fun, sep="")
        source(source_fun) # sourcing fun into R
        
        
        seq_regs <- as.character(segmented_file[,4])
        
        
        
        #######Plots#######  
        ######CNA PROFILE HEATMAP#####
        output$cna_profile_heatmap_MAIN <- renderUI({
          
          HTML(paste("<br><h4><b>Copy number genome profiles</b></h4>(samples ordered by <i><b>", re_order_var, "</b></i>)", sep=""))
        })  
        ## Sample cna profile heatmap ##
        data <- data.means
        colnames(data) <- sample_names
        rownames(data) <- seg_regions_names
        
        mat_vars <- mat_variables
        track_vars <- var_tracks
        gain.cutoff <- re_gain_th
        loss.cutoff <- re_loss_th
        
        # Preparing p_annotation (Annotation tracks):
        fun_name <- "Ht.Annotation.RegionProfile" # function name
        my_fun <- paste(fun_name, ".R", sep="") # file function name
        source_fun <- paste(dir_funs, "/", my_fun, sep="")
        source(source_fun) # sourcing fun into R
        
        p_annotation <- Ht.Annotation.RegionProfile(data, mat_vars, variables_info_mat, track_vars)
        ############################################
        
        # Preparing p_annotation2 (to join with p_chroms):
        fun_name <- "Ht.Annotation.2.RegionProfile" # function name
        my_fun <- paste(fun_name, ".R", sep="") # file function name
        source_fun <- paste(dir_funs, "/", my_fun, sep="")
        source(source_fun) # sourcing fun into R
        
        p_annotation2 <- Ht.Annotation.2.RegionProfile(data, mat_vars, variables_info_mat, track_vars)
        ############################################
        
        # Preparing p_chroms:
        fun_name <- "Ht.Chroms.RegionProfile" # function name
        my_fun <- paste(fun_name, ".R", sep="") # file function name
        source_fun <- paste(dir_funs, "/", my_fun, sep="")
        source(source_fun) # sourcing fun into R
        
        p_chroms <- Ht.Chroms.RegionProfile(segmented_file)
        ############################################
        
        # Preparing p_cn (nude Region Profile heatmap):
        fun_name <- "Ht.RegionProfile" # function name
        my_fun <- paste(fun_name, ".R", sep="") # file function name
        source_fun <- paste(dir_funs, "/", my_fun, sep="")
        source(source_fun) # sourcing fun into R
        
        p_cna_profile <- Ht.RegionProfile(data, mat_vars, segmented_file, variables_info_mat, track_vars, gain.cutoff, loss.cutoff)
        ############################################
        
        # Merging final heatmap:
        h1 <- (0.3*length(var_tracks))/6
        h2 <- 1-h1
        seq_heights <- c(h1, h2)
        seq_widths <- c(0.05, 0.95)
        
        pp <- subplot(p_annotation, p_cna_profile, nrows=2, shareX=T, heights=seq_heights)
        pp2 <- subplot(p_annotation2, p_chroms, nrows=2, shareX=T, heights=seq_heights)
        p_cn <- subplot(pp2, pp, shareY=F, widths=seq_widths)
        ############################################
        
        
        
        ## Sample cna profile heatmap ##
        
        
        output$cna_profile_heatmap <- renderPlotly({
          p_cn
          
        })
        # download buttons:
        output$dw_buttons_cna_profile_htmp <- renderUI({
          div(style="text-align:right;",
              downloadButton("TSV_cna_profile_htmp", "TSV"),
              downloadButton("PNG_cna_profile_htmp", "PNG")
          )
        })
        
        
        ######CNA PROFILE dendrogram euclidian distance#####
        output$cna_profile_dendro_MAIN <- renderUI({
          HTML("<h4 style='text-align:center'><b>Hierarchical clustering for CN region profile</b></h4>")
        })
        
        output$cna_prof_dendro_plot_button <- renderUI({
          div(
            ##Button run cin
            div(style="text-align:center;padding-top:2em;", actionBttn("button_cnaProfDendro", "Plot!", style="simple", size="sm", color="primary")),
            busyIndicator(text="", wait=200)
          )
        })
        
        observeEvent(input$button_cnaProfDendro,{
          
          ####the plot####
          ## Sample cna profile DENDROGRAM ##
          
          
          # output$dendro_type_box <- renderUI({
          #   div(
          #     selectizeInput("dendro_type", "Dendrogram type:", choices = c("both", "column", "row"), selected= "row", width="200px")
          #   )
          # })
          
          # Preparing p_dendro:
          fun_name <- "dendro.Ht.RegionProfile" # function name
          my_fun <- paste(fun_name, ".R", sep="") # file function name
          source_fun <- paste(dir_funs, "/", my_fun, sep="")
          source(source_fun) # sourcing fun into R
          
          data <- data.means
          colnames(data) <- sample_names
          rownames(data) <- seg_regions_names
          
          mat_vars <- mat_variables
          track_vars <- var_tracks
          gain.cutoff <- gain_th
          loss.cutoff <- loss_th
          
          # p_dendro <<- reactive({
          
          dendro <- "both"
          # dendro <- as.character(input$dendro_type)
          
          p_dendro <- dendro.Ht.RegionProfile(data, dendro, mat_vars, segmented_file, variables_info_mat, track_vars, gain.cutoff, loss.cutoff)
          
          # })
          
          output$prova_1 <- renderUI({
            HTML("dimensions mat ", dim(mat))
          })
          output$prova_2 <- renderUI({
            HTML("dimensions annot ", dim(annot))
          })
          output$prova_3 <- renderUI({
            HTML("sample variable ", sample_variable)
          })
          output$prova_4 <- renderUI({
            HTML("dimensions mat_variables ", length(mat_variables[,sample_variable]))
          })
          
          output$cna_profile_dendrogram <- renderPlotly({
            p_dendro
            
          })
          
          ## Sample cna profile DENDROGRAM ##
          
          # download buttons:
          output$dw_buttons_dendro_profile <- renderUI({
            div(style="text-align:right;",
                downloadButton("PNG_dendro_profile", "PNG")
            )
          })
          
        })
        
        
        
        #####Genes in selected region (below cna profile)####
        
        output$select_reg_genes <- renderUI({
          div(
            selectizeInput("reg_selected_in_cn_tab", "", choices = seg_regions_names, options = list(maxOptions=3000), selected= seg_regions_names[1])
          )
        })
        
        reg_selected <- reactive({
          term <- as.character(input$reg_selected_in_cn_tab)
          term
        })
        output$tab_genes_MAIN <- renderUI({
          HTML("<br><h4>Genes in region: <i style=color:orange>(by ", segmentation_name, ")</i>")
        })
        
        
        tab_genes <- reactive({
          term <- as.character(input$reg_selected_in_cn_tab)
          term
          
          row_region <- segmented_file[which(segmented_file[,4]==term),]
          
          fun_name <- "refGene.selection" # function name
          my_fun <- paste(fun_name, ".R", sep="")
          source_fun <- paste(dir_funs, "/", my_fun, sep="")
          source(source_fun)
          # sourcing fun into R
          
          refGene_mat <- refGene.selection(genome_build)
          
          CHROM <- as.character(row_region[,1])
          START <- as.numeric(as.character(row_region[,2]))
          END <- as.numeric(as.character(row_region[,3]))
          
          rows_in_refGene <- which(refGene_mat[,"chrom"]==CHROM & as.numeric(as.character(refGene_mat[,"start"]))>=START & as.numeric(as.character(refGene_mat[,"end"]))<=END)
          
          
          
          genes_in_reg <- refGene_mat[rows_in_refGene,c("symbol_name", "chrom", "strand", "start", "end")]
          uni_g <- unique(as.character(genes_in_reg[,1]))
          each_gene <- as.vector(do.call(cbind, mclapply(uni_g, function(x,y){which(y==x)[1]}, y=as.character(genes_in_reg[,1]))))
          
          genes_in_reg_2 <- genes_in_reg[each_gene,]
          genes_in_reg_2
        })
        
        output$tab_genes_reg <- renderDataTable({
          return(tab_genes())
        })
        # download buttons:
        output$dw_tab_genes_bttn <- renderUI({
          div(style="text-align:left;",
              downloadButton("TSV_tab_genes_in_cn_profile", "TSV")
          )
        })
        
        #####Unsupervised CORRELATION SAMPLES HEATMAP####
        
        output$unsupervised_corr_heatmap_MAIN <- renderUI({
          HTML(paste("<br><h4><b>Samples correlation heatmap</b></h4>(samples ordered by <i><b>", order_var_name, "</b></i>)<br>(correlation method = <i><b>", re_method_cor, "</b></i>)", sep=""))
        })  
        ####the plot####
        
        # Preparing p_corr_samples:
        fun_name <- "CorrSamples.Ht.RegionProfile" # function name
        my_fun <- paste(fun_name, ".R", sep="") # file function name
        source_fun <- paste(dir_funs, "/", my_fun, sep="")
        source(source_fun) # sourcing fun into R
        
        data <- data.means
        colnames(data) <- sample_names
        rownames(data) <- seg_regions_names
        
        corr.method <- re_method_cor
        
        p_corr_samples <- CorrSamples.Ht.RegionProfile(data, corr.method)
        ############################
        
        # Merging final p_cor:
        h1 <- (0.35*length(var_tracks))/6
        h2 <- 1-h1
        seq_heights <- c(h1, h2)
        
        p_corr <- subplot(p_annotation, p_corr_samples, nrows=2, shareX=T, heights=seq_heights)
        ############################
        
        
        output$unsupervised_corr_sample_heatmap <- renderPlotly({
          ##CORRELATION HEATMAP
          p_corr
          
          ## Unsupervised sample correlation heatmap ##
          
        })
        # download buttons:
        output$dw_buttons_unsup_corr_htmp <- renderUI({
          div(style="text-align:right;",
              downloadButton("TSV_unsup_corr_htmp", "TSV"),
              downloadButton("PNG_unsup_corr_htmp", "PNG")
          )
        })
        
        ##Annotation tracks buttons dw:
        output$dw_annotation_tracks_tab1 <- renderUI({
          div(style="text-align:left;",
              downloadButton("TSV_annotation_tracks_tab1", "Annotation tracks TSV")
          )
        })
        output$dw_annotation_tracks_tab2 <- renderUI({
          div(style="text-align:left;",
              downloadButton("TSV_annotation_tracks_tab2", "Annotation tracks TSV")
          )
        })
        output$dw_annotation_tracks_tab3 <- renderUI({
          div(style="text-align:left;",
              downloadButton("TSV_annotation_tracks_tab3", "Annotation tracks TSV")
          )
        })
        
        
        
        
        
        
        
        ######CORRELATION dendrogram euclidian distance#####
        output$cor_samples_dendro_MAIN <- renderUI({
          HTML(paste("<b><h4><b>Hierarchical clustering for samples correlation</b></h4></b>(correlation method = <i><b>", re_method_cor, "</b></i>)", sep=""))
        })  
        
        output$dendro_type_box_cor_plot_button <- renderUI({
          div(
            ##Button run cin
            div(style="text-align:center;padding-top:2em;", actionBttn("button_dendroCor", "Plot!", style="simple", size="sm", color="primary")),
            busyIndicator(text="", wait=200)
          )
        })
        
        
        observeEvent(input$button_dendroCor,{
          ####the plot####
          ## Sample cna profile DENDROGRAM ##
          
          fun_name <- "dendro.Ht.CorrSamples.RegionProfile" # function name
          my_fun <- paste(fun_name, ".R", sep="") # file function name
          source_fun <- paste(dir_funs, "/", my_fun, sep="")
          source(source_fun) # sourcing fun into R
          
          data <- data.means
          colnames(data) <- sample_names
          rownames(data) <- seg_regions_names
          
          corr.method <- re_method_cor
          mat_vars <- mat_variables
          track_vars <- var_tracks
          gain.cutoff <- gain_th
          loss.cutoff <- loss_th
          
          p_dendro_cor <<- dendro.Ht.CorrSamples.RegionProfile(data, corr.method, mat_vars, variables_info_mat, track_vars, gain.cutoff, loss.cutoff)
          
          
          output$cor_samples_dendrogram <- renderPlotly({
            p_dendro_cor
            
          })
          # download buttons:
          output$dw_buttons_dendro_cor <- renderUI({
            div(style="text-align:right;",
                downloadButton("PNG_dendro_cor", "PNG")
            )
          })
        })
        
        
        
        ## Sample cna profile DENDROGRAM ##
        ####CNA-regions EVENTS####
        colnames(data.means) <- sample_names
        
        fun_name <- "fun.grep" # function name
        my_fun <- paste(fun_name, ".R", sep="") # file function name
        source_fun <- paste(dir_funs, "/", my_fun, sep="")
        source(source_fun) # sourcing fun into R
        
        
        fun_name <- "extract.stats" # function name
        my_fun <- paste(fun_name, ".R", sep="") # file function name
        source_fun <- paste(dir_funs, "/", my_fun, sep="")
        source(source_fun) # sourcing fun into R
        
        
        seq_regs <- as.character(segmented_file[,4])
        
        ## Statitstics CNA counts by samples
        events_order <- c("normal", "loss", "gain")
        counts_0 <- as.matrix(apply(data.events, 2, fun_grep, z=events_order))
        # counts_0_freq <- round((counts_0/nrow(segmented_file) ) *100,2) # counts as frequency percentage!
        # counts <- as.data.frame(cbind(sample_names,t(counts_0_freq))) 
        counts <- as.data.frame(cbind(sample_names,t(counts_0)))
        colnames(counts)[1] <- "samples"
        
        counts_by_samples <- counts
        rownames(counts_by_samples) <- sample_names
        
        ans_gain_byS <- extract.stats(counts_by_samples, "gain", "samples")
        ans_loss_byS <- extract.stats(counts_by_samples, "loss", "samples")
        
        ## Statitstics CNA counts by region
        events_order <- c("normal", "loss", "gain")
        counts_0 <- as.matrix(apply(data.events, 1, fun_grep, z=events_order))
        counts_0_freq <- round((counts_0/length(sample_names) ) *100,2)
        counts <- as.data.frame(cbind(seq_regs,t(counts_0_freq)))
        colnames(counts)[1] <- "regions"
        
        counts_by_regions <- counts
        rownames(counts_by_regions) <- seq_regs
        
        ans_gain_byR <- extract.stats(counts_by_regions, "gain", "regions")
        ans_loss_byR <- extract.stats(counts_by_regions, "loss", "regions")
        
        ## Statitstics CNA region means by sample
        
        xx <- cbind(data.means, seq_regs)
        
        ans <- mclapply(sample_names, function(x){extract.stats(xx, x, "seq_regs")})
        names(ans) <- sample_names
        final_ans <- ans
        
        min_values_samples <- as.vector(do.call(rbind, mclapply(final_ans, function(x) as.numeric(as.character(x[["min"]])))))
        top3_MIN_samples <- sample_names[order(min_values_samples)][1:3]
        top3_MIN_samples_values <- min_values_samples[order(min_values_samples)][1:3]
        
        max_values_samples <- as.vector(do.call(rbind, mclapply(final_ans, function(x) as.numeric(as.character(x[["max"]])))))
        top3_MAX_samples <- rev(sample_names[order(max_values_samples)])[1:3]
        top3_MAX_samples_values <- rev(max_values_samples[order(max_values_samples)])[1:3]
        
        ## Statitstics CNA region means by region
        
        xx <- cbind(t(data.means),sample_names)
        rownames(xx) <- sample_names
        
        seq_regs <- as.character(segmented_file[,4])
        
        ans <- mclapply(seq_regs, function(x){extract.stats(xx, x, "sample_names")})
        names(ans) <- seq_regs
        final_ans <- ans
        
        min_values_regions <- as.vector(do.call(rbind, mclapply(final_ans, function(x) as.numeric(as.character(x[["min"]])))))
        top3_MIN_region <- seq_regs[order(min_values_regions)][1:3]
        top3_MIN_region_values <- min_values_regions[order(min_values_regions)][1:3]
        
        max_values_regions <- as.vector(do.call(rbind, mclapply(final_ans, function(x) as.numeric(as.character(x[["max"]])))))
        top3_MAX_region <- rev(seq_regs[order(max_values_regions)])[1:3]
        top3_MAX_region_values <- rev(max_values_regions[order(max_values_regions)])[1:3]
        
        #####Bar plot by samples#######
        
        
        ############################################
        
        
        output$cna_events_MAIN <- renderUI({
          HTML(paste("<br><h4><b>CNAs region counts by samples</b></h4>(samples ordered by <i><b>", re_order_var, "</b></i>)", sep=""))
        })
        # the plot:
        # Preparing 'p_events' by samples:
        fun_name <- "freq.barplot" # function name
        my_fun <- paste(fun_name, ".R", sep="")
        source_fun <- paste(dir_funs, "/", my_fun, sep="")
        source(source_fun)
        # sourcing fun into R
        
        colors <- c("blue", "gainsboro", "red")
        data <- as.data.frame(counts_by_samples[,c("loss", "normal", "gain")])
        
        p_bar <- freq.barplot(data,colors)
        
        h1 <- (0.35*length(var_tracks))/6
        h2 <- 1-h1
        seq_heights <- c(h1, h2)
        
        p_events <- subplot(p_annotation, p_bar, nrows=2, shareX=T, heights=seq_heights, titleY=T)
        
        output$cna_events <- renderPlotly({
          p_events
        })
        #download buttons:
        output$dw_buttons_cna_events_barplot <- renderUI({
          div(style="text-align:right;",
              downloadButton("TSV_cna_events_barplot", "TSV"),
              downloadButton("PNG_cna_events_barplot", "PNG")
          )
        })
        
        #####Bar plot by regions#######
        
        # Preparing 'p_bar_by_regs' by regions:
        fun_name <- "freq.barplot.GainLoss" # function name
        my_fun <- paste(fun_name, ".R", sep="")
        source_fun <- paste(dir_funs, "/", my_fun, sep="")
        source(source_fun)
        # sourcing fun into R
        
        colors <- c("blue", "gainsboro", "red")
        data <- as.data.frame(counts_by_regions[,c("loss", "normal", "gain")])
        
        p_bar_by_regs <- freq.barplot.GainLoss(data, colors, segmented_file)
        ############################################
        
        
        output$cna_events_by_regs_MAIN <- renderUI({
          HTML(paste("<br><h4><b>CNA region frequencies</b></h4>", sep=""))
        })
        # the plot:
        output$cna_events_by_regs <- renderPlotly({
          
          ## CNA events plot summary ##
          p_bar_by_regs
          
          ## CNA events plot summary ##
          #####################################
        })
        #download buttons:
        output$dw_buttons_cna_events_by_regs <- renderUI({
          div(style="text-align:right;",
              downloadButton("TSV_cna_events_by_regs", "TSV"),
              downloadButton("PNG_cna_events_by_regs", "PNG")
          )
        })
        
        #########################################################
        #########################################################
        #######################DOWNLOAD HANDLER:####################
        ########TSV_annotation_tracks_tab1#####
        output$TSV_annotation_tracks_tab1 <- downloadHandler(
          filename = paste("Annotation_tracks_by_", re_order_var, "_", Sys.time(), ".tsv", sep=""),
          content = function(file) {
            write.table(mat_variables, file, sep="\t", col.names=NA, row.names=T, quote=F, dec=".")
          }
        )
        ########TSV_annotation_tracks_tab2#####
        output$TSV_annotation_tracks_tab2 <- downloadHandler(
          filename = paste("Annotation_tracks_by_", re_order_var, "_", Sys.time(), ".tsv", sep=""),
          content = function(file) {
            write.table(mat_variables, file, sep="\t", col.names=NA, row.names=T, quote=F, dec=".")
          }
        )
        ########TSV_annotation_tracks_tab3#####
        output$TSV_annotation_tracks_tab3 <- downloadHandler(
          filename = paste("Annotation_tracks_by_", re_order_var, "_", Sys.time(), ".tsv", sep=""),
          content = function(file) {
            write.table(mat_variables, file, sep="\t", col.names=NA, row.names=T, quote=F, dec=".")
          }
        )
        ########TSV_unsup_corr_htmp#####
        output$TSV_unsup_corr_htmp <- downloadHandler(
          filename = paste("Unsup_corr_mat_", segmentation_name, "_by_", re_order_var, "_", Sys.time(), ".tsv", sep=""),
          content = function(file) {
            mat <- data.means
            colnames(mat) <- sample_names
            cor_samples <- cor(mat, method=re_method_cor) # 'method_cor' is selected by user (pearson, kendall, spearman)
            cor_samples_mat <- round(cor_samples, 3) # sample correlations matrix
            
            write.table(cor_samples_mat, file, sep="\t", col.names=NA, row.names=T, quote=F, dec=".")
          }
        )
        ########PNG_unsup_corr_htmp#####
        output$PNG_unsup_corr_htmp <- downloadHandler(
          filename = paste("Unsup_corr_heatmap_", segmentation_name, "_", re_method_cor, "_", Sys.time(), ".png", sep=""),
          content = function(file) {
            p_corr$x$layout$font$size <- 70
            p_corr$x$layout$margin$t <- 400
            p_corr$x$layout$margin$b <- 400
            p_corr$x$layout$margin$l <- 700
            p_corr$x$layout$margin$r <- 700
            export(p_corr, file = file, vwidth = 7000, vheight = 5000, expand=1)
            
          }
        )
        ########TSV_cna_profile_htmp#####
        colnames(data.means) <- sample_names
        rownames(data.means) <- seg_regions_names
        
        
        output$TSV_cna_profile_htmp <- downloadHandler(
          filename = paste("cna_profile_", segmentation_name,  "_by_", re_order_var, Sys.time(), ".tsv", sep=""),
          content = function(file) {
            write.table(data.means, file, sep="\t", col.names=NA, row.names=T, quote=F, dec=".")
          }
        )
        ########PNG_cna_profile_htmp#####
        output$PNG_cna_profile_htmp <- downloadHandler(
          filename = paste("cna_profile_", segmentation_name, "_by_", re_order_var, "_", Sys.time(), ".png", sep=""),
          content = function(file) {
            p_cn$x$layout$font$size <- 70
            p_cn$x$layout$margin$t <- 700
            p_cn$x$layout$margin$b <- 1100
            p_cn$x$layout$margin$l <- 1000
            p_cn$x$layout$margin$r <- 700
            export(p_cn,file = file, vwidth = 7000, vheight = 5000, expand=1)
            
          }
        )
        
        
        
        
        ########TSV_cna_events_barplot#####
        output$TSV_cna_events_barplot <- downloadHandler(
          filename = paste("cna_events_in_", segmentation_name, "_by_", re_order_var, "_", Sys.time(), ".tsv", sep=""),
          content = function(file) {
            write.table(counts_by_samples, file, sep="\t", col.names=T, row.names=F, quote=F, dec=".")
          }
        )
        ########PNG_cna_events_barplot#####
        output$PNG_cna_events_barplot <- downloadHandler(
          filename = paste("cna_events_in_", segmentation_name, "_by_", re_order_var, "_", Sys.time(), ".png", sep=""),
          content = function(file) {
            p_events$x$layout$font$size <- 70
            p_events$x$layout$margin$t <- 500
            p_events$x$layout$margin$b <- 500
            p_events$x$layout$margin$l <- 500
            p_events$x$layout$margin$r <- 500
            export(p_events, file = file, vwidth = 7000, vheight = 5000, expand=1)
            
          }
        )
        ########TSV_cna_events_barplot#####
        output$TSV_cna_events_by_regs <- downloadHandler(
          filename = paste("cna_frequencies_in_", segmentation_name, "_by_", re_order_var, "_", Sys.time(), ".tsv", sep=""),
          content = function(file) {
            write.table(counts_by_regions, file, sep="\t", col.names=T, row.names=F, quote=F, dec=".")
          }
        )
        ########PNG_cna_events_barplot#####
        output$PNG_cna_events_by_regs <- downloadHandler(
          filename = paste("region_frequencies_", segmentation_name, "_", Sys.time(), ".png", sep=""),
          content = function(file) {
            if ( nrow(segmented_file)>88 ){
              p_bar_by_regs$x$layout$bargap <- 0
            }
            p_bar_by_regs$x$layout$font$size <- 80
            p_bar_by_regs$x$layout$xaxis$tickfont$size <- 100
            p_bar_by_regs$x$layout$margin$t <- 300
            p_bar_by_regs$x$layout$margin$b <- 700
            p_bar_by_regs$x$layout$margin$l <- 500
            p_bar_by_regs$x$layout$margin$r <- 500
            export(p_bar_by_regs, file = file, vwidth = 7000, vheight = 5000, expand=1)
          }
        )
        ########TSV_tab_genes_in_cn_profile#####
        output$TSV_tab_genes_in_cn_profile <- downloadHandler(
          filename = reactive({
            term <- as.character(input$reg_selected_in_cn_tab)
            paste("Genes_in_region_", as.character(term), "_", segmentation_name, "_", Sys.time(), ".tsv", sep="")
          }),
          content = function(file) {
            write.table(as.data.frame(tab_genes()), file, sep="\t", col.names=T, row.names=F, quote=F, dec=".")
          }
        )
        
        
        ########################## #
        
        
        
        
        
        ########PNG_dendro_cor#####
        output$PNG_dendro_cor <- downloadHandler(
          filename = paste("Hierarchical_clust_corr_samples_", Sys.time(), ".png", sep=""),
          content = function(file) {
            
            
            export(p_dendro_cor, file = file, vwidth = 1000, vheight = 750, expand=1)
          }
        )
        ########PNG_dendro_profile#####
        output$PNG_dendro_profile <- downloadHandler(
          filename = paste("Hierarchical_clust_profile_", Sys.time(), ".png", sep=""),
          content = function(file) {
            data <- data.means
            colnames(data) <- sample_names
            rownames(data) <- seg_regions_names
            
            mat_vars <- mat_variables
            track_vars <- var_tracks
            gain.cutoff <- gain_th
            loss.cutoff <- loss_th
            
            # p_dendro <<- reactive({
            
            dendro <- "both"
            # dendro <- as.character(input$dendro_type)
            
            p_dendro <-  dendro.Ht.RegionProfile(data, dendro, mat_vars, segmented_file, variables_info_mat, track_vars, gain.cutoff, loss.cutoff, fontsize_row=10)
            
            # })
            
            export(p_dendro, file = file, vwidth = 1000, vheight = 750, expand=1)
          }
        )
        
        
        
        
        
      } # else for (length(var_tracks)<=6)
      
      
      #########################################################
      #########################################################
      #########################################################
    })
    
    ##########WHEN BUTTON RUN LIMMA###########
    observeEvent(input$button_run_limma, {
      
      
      ###############TESTS t-test & fisher part######################
      
      ## variables control ##
      
      p_val_type <- as.character(input$pval_type)
      
      
      group_variable <- input$group_var
      class_var <- as.character(variables_info_mat[which(variables_info_mat==group_variable),"class_var"])
      if (class_var=="numeric") {
        values1 <- as.numeric(as.character(input$slider_group_var))
        min1 <- min(values1, na.rm = T)
        max1 <- max(values1, na.rm = T)
        rows_in_annotation1 <- which(as.numeric(as.character(mat_variables[,group_variable]))>=min1 & as.numeric(as.character(mat_variables[,group_variable]))<=max1)
        
        values2 <- as.numeric(as.character(input$slider_group_var_2))
        min2 <- min(values2, na.rm = T)
        max2 <- max(values2, na.rm = T)
        rows_in_annotation2 <- which(as.numeric(as.character(mat_variables[,group_variable]))>=min2 & as.numeric(as.character(mat_variables[,group_variable]))<=max2)
        
        rows_in_annotation <- unique(c(rows_in_annotation1, rows_in_annotation2))
        col_group_var <- sort(as.numeric(as.character(mat_variables[rows_in_annotation,group_variable])))
        real_groups <- unique(col_group_var)
      } 
      if (class_var=="categoric"){
        real_groups <- input$check_group_var
        rows_in_annotation <- which(mat_variables[,group_variable]%in%real_groups)
        col_group_var <- as.character(mat_variables[rows_in_annotation,group_variable])
      }
      n_groups <- length(real_groups)
      
      min_pval <- min(as.numeric(as.character(input$slider_pval_limma)))
      max_pval <- max(as.numeric(as.character(input$slider_pval_limma)))
      
      data_0 <- data.means[,rows_in_annotation]
      rownames(data_0) <- seg_regions_names
      
      if (n_groups<=1) {
        output$text_limma <- renderUI(div(style="color:red","You have to select MORE groups!!"))
        
        output$limma_heatmap_MAIN <- renderUI({
          HTML("")
        })
        output$limma_heatmap <- renderPlotly({
          
          plotly_empty()
          
        })
        output$reg_boxplot_MAIN <- renderUI({
          HTML("")
        })
        output$reg_boxplot <- renderPlotly({
          
          plotly_empty()
          
        })
        
        output$dw_buttons_limma_heatmap <- renderUI({
        })
        output$dw_buttons_reg_boxplot <- renderUI({
        })
        
        output$not_differences_2 <- output$not_differences <- renderUI(div(""))
        
        
        ###
        
        output$fisher_heatmap_MAIN <- renderUI({
          HTML("")
        })
        output$fisher_heatmap <- renderPlotly({
          
          plotly_empty()
          
        })
        output$reg_boxplot_MAIN_fisher <- renderUI({
          HTML("")
        })
        output$reg_boxplot_fisher <- renderPlotly({
          
          plotly_empty()
          
        })
        
        output$dw_buttons_fisher_heatmap <- renderUI({
        })
        output$dw_buttons_reg_boxplot_fisher <- renderUI({
        })
        
        output$not_differences_2_fisher <- output$not_differences_fisher <- renderUI(div(""))
        
      } else if (n_groups>35) {
        output$text_limma <- renderUI(div(style="color:red","Too many groups!! Please, reduce selection."))
        
        output$limma_heatmap_MAIN <- renderUI({
          HTML("")
        })
        output$limma_heatmap <- renderPlotly({
          
          plotly_empty()
          
        })
        output$reg_boxplot_MAIN <- renderUI({
          HTML("")
        })
        output$reg_boxplot <- renderPlotly({
          
          plotly_empty()
          
        })
        
        output$dw_buttons_limma_heatmap <- renderUI({
        })
        output$dw_buttons_reg_boxplot <- renderUI({
        })
        
        output$not_differences_2 <- output$not_differences <- renderUI(div(""))
        
        
        ###
        
        output$fisher_heatmap_MAIN <- renderUI({
          HTML("")
        })
        output$fisher_heatmap <- renderPlotly({
          
          plotly_empty()
          
        })
        output$reg_boxplot_MAIN_fisher <- renderUI({
          HTML("")
        })
        output$reg_boxplot_fisher <- renderPlotly({
          
          plotly_empty()
          
        })
        
        output$dw_buttons_fisher_heatmap <- renderUI({
        })
        output$dw_buttons_reg_boxplot_fisher <- renderUI({
        })
        
        output$not_differences_2_fisher <- output$not_differences_fisher <- renderUI(div(""))
        
      } else if (length(col_group_var)==length(real_groups)) {
        output$text_limma <- renderUI(div(style="color:red","Variable cannot define groups (non replicative samples...). Choose another one!!"))
        
        output$limma_heatmap_MAIN <- renderUI({
          HTML("")
        })
        output$limma_heatmap <- renderPlotly({
          
          plotly_empty()
          
        })
        output$reg_boxplot_MAIN <- renderUI({
          HTML("")
        })
        output$reg_boxplot <- renderPlotly({
          
          plotly_empty()
          
        })
        
        output$dw_buttons_limma_heatmap <- renderUI({
        })
        output$dw_buttons_reg_boxplot <- renderUI({
        })
        
        output$not_differences_2 <- output$not_differences <- renderUI(div(""))
        
        
        ###
        
        output$fisher_heatmap_MAIN <- renderUI({
          HTML("")
        })
        output$fisher_heatmap <- renderPlotly({
          
          plotly_empty()
          
        })
        output$reg_boxplot_MAIN_fisher <- renderUI({
          HTML("")
        })
        output$reg_boxplot_fisher <- renderPlotly({
          
          plotly_empty()
          
        })
        
        output$dw_buttons_fisher_heatmap <- renderUI({
        })
        output$dw_buttons_reg_boxplot_fisher <- renderUI({
        })
        
        output$not_differences_2_fisher <- output$not_differences_fisher <- renderUI(div(""))
        
      } else if (n_groups >1 ) {
        output$text_limma <- renderUI(paste("You have selected", n_groups, "groups!", sep=" "))
        
        ############ Boxplot by region############
        
        
        reg_to_plot <- reactive({
          term <- as.character(input$reg_selected)
          term
        })
        
        
        
        ### T-STUDENT PART (by limma) ####
        
        groups <- paste("Group_", 1:length(real_groups), sep="")
        
        
        mm <- model.matrix(~0+factor(col_group_var))
        colnames(mm) <- groups
        
        e_vector <- rep(NA, length=1)
        for (i in 1:length(real_groups)){
          for (j in (i:length(real_groups))[-(i:length(real_groups)==i)]){
            contr <- paste(real_groups[j],real_groups[i],sep="-")
            e_vector <- c(e_vector, contr)
          }
        }
        
        real_CONTRs <- e_vector[-1]
        # print(summary(mat_variables.0))
        e_vector <- rep(NA, length=1)
        for (i in 1:length(groups)){
          for (j in (i:length(groups))[-(i:length(groups)==i)]){
            contr <- paste(groups[j],groups[i],sep="-")
            e_vector <- c(e_vector, contr)
          }
        }
        
        CONTRs <- e_vector[-1]
        
        
        data_0 <- data.means[,rows_in_annotation]
        rownames(data_0) <- seg_regions_names
        not_sd <- which(apply(data_0, 1, sd)==0)
        #data_1 <- data_0[-not_sd,]
        data_1 <- data_0
        fit <- lmFit(data_1, mm)
        contrast.matrix <- makeContrasts(contrasts=CONTRs, levels=mm)
        fit2 <- contrasts.fit(fit, contrast.matrix)
        fit2.1 <- eBayes(fit2)
        
        
        if (p_val_type == "p-value") {
          pval_term <- "p-val"
          col_values <- 4
        } else if (p_val_type == "Adj.p-value") {
          pval_term <- "Adj.p-val"
          col_values <- 5
        }
        
        
        e_list <- list()
        for (i in 1:length(CONTRs)){
          tab_0 <- topTable(fit2.1, coef=i, number=nrow(data_1))
          e_vector <- rep(NA, length(seg_regions_names))
          for (j in 1:length(seg_regions_names)){
            reg <- seg_regions_names[j]
            row <- which(rownames(tab_0)==reg)
            if(length(row)==1){
              p_value <- tab_0[row,col_values]
            } else if (length(row)==0) {
              p_value <- 1
            }
            e_vector[j] <- p_value
          }
          
          e_list[[paste(CONTRs[i])]] <- e_vector
          
          #	NAME <- paste(out_dir, "lmTAB_", name_SEG_FILE, "_", CONTRs[i], ".txt", sep="")
          #	write.table(tab, NAME, sep="\t", col.names=NA, row.names=T, quote=F)
        }
        
        limma_df <- as.data.frame(as.matrix(do.call(cbind, e_list)))
        colnames(limma_df) <- CONTRs
        rownames(limma_df) <- seg_regions_names
        
        limma_df_2 <- limma_df
        colnames(limma_df_2) <- real_CONTRs
        
        #	NAME <- paste(out_dir, "more_diff_regions_adj.Pval_",  name_SEG_FILE, ".txt", sep="")
        #	write.table(limma_df, NAME, sep="\t", col.names=NA, row.names=T, quote=F)
        
        
        annot_limma <- as.data.frame(as.matrix(cbind(CONTRs, rep(NA, length(CONTRs)))))
        annot_limma[,1] <- do.call(rbind,strsplit(CONTRs, "-"))[,2]
        annot_limma[,2] <- do.call(rbind,strsplit(CONTRs, "-"))[,1]
        rownames(annot_limma) <- CONTRs
        colnames(annot_limma) <- c("Group.1","Group.2")
        for (i in 1:nrow(annot_limma)){
          for (j in 1:ncol(annot_limma)){
            cell <- annot_limma[i,j]
            annot_limma[i,j] <- as.character(real_groups[which(groups==cell)])
            
          }
        }
        
        ## Limma analysis ##
        
        
        
        
        ####### OUTPUTS Heatmap regions########
        # LIMMA COMPARISIONS HEATMAP:
        
        output$limma_heatmap_MAIN <- renderUI({
          HTML(paste("<br><h4>More differentiated regions between groups defined by <i><b>",
                     group_variable, "</b></i>", sep=""))
          # withTags(
          #   hr(),
          #   h4(paste("More differentiated regions between groups defined by ", tag$b(group_variable), sep="")),
          #   paste("(adjusted p-value)", sep="")
          #   
          # )
        })
        # the plot:
        n.colors <- length(seq(by=0.01, from=-2, to=2)) # centrat a 0
        #	n.colors <- round((MAX-MIN)+1)+100 # centrat a la mediana de la sèrie entre MAX i MIN
        colfunc <- colorRampPalette(c("red", "yellow", "white"))(n.colors)
        
        mat2 <- as.matrix(limma_df)
        # colnames(mat2) <- CONTRs
        #rownames(mat2) <- paste(segmented_file[,1], segmented_file[,4])
        for (i in 1:nrow(mat2)) {
          for (j in 1:ncol(mat2)) {
            #if (mat2[i,j]<(0.05)) mat2[i,j]<- (0.05)
            if (mat2[i,j]>(max_pval)) mat2[i,j]<- (max_pval)
            if (mat2[i,j]<(min_pval)) mat2[i,j]<- (min_pval)
            
          }
        }
        
        if (length(which(apply(mat2, 2, sd)!=0))==0){
          #output$limma_heatmap_MAIN <- renderUI({HTML("")})
          #output$text_limma <- renderUI(div(style="color:red","Try to expand p-value range!!"))
          output$not_differences_2 <- output$not_differences <- renderUI(HTML(paste("<p style='color:red'>Sorry, but there are no descriptive regions between '", group_variable, "' groups<br>Try to expand p-value range!!</p>")))
          output$limma_heatmap <- renderPlotly({
            
            plotly_empty()
            
          })
          # output$dw_buttons_limma_heatmap <- renderUI({
          #   div(style="text-align:right;",
          #       actionButton("empty_1", "TSV"),
          #       actionButton("empty_2", "PNG")
          #   )
          # })
          output$dw_buttons_reg_boxplot <- renderUI({
            div(style="text-align:right;",
                actionButton("empty_3", "TSV"),
                actionButton("empty_4", "PNG")
            )
          })
        } else if (length(which(apply(mat2, 2, sd)!=0))>0) {
          output$not_differences_2 <- output$not_differences <- renderUI(div(""))
          
          ####the plot####
          
          
          m <- as.matrix(t(annot_limma))
          m1 <- mapply(m, FUN=as.factor)
          data.groups <- matrix(data=as.numeric(m1), ncol=ncol(m), nrow=nrow(m))
          colnames(data.groups) <- real_CONTRs
          rownames(data.groups) <- colnames(annot_limma)
          
          y <- rownames(data.groups)
          palette <- "Spectral"
          palette <- as.character(variables_info_mat[which(variables_info_mat[,"name_var"]==group_variable),"color_palette"])
          
          showscale <- TRUE
          x <- paste("'", as.character(colnames(data.groups)), "'", sep="")
          
          text <- t(m)
          
          mat_test <- t(data.groups)
          tickvals <- unique(as.numeric(m1))
          range <- tickvals[length(tickvals)]-tickvals[1]
          tickvals[1] <- tickvals[1] + (range*0.1)
          tickvals[length(tickvals)] <- tickvals[length(tickvals)] - (range*0.1)
          ticktext <- levels(m1)
          
          
          p_groups <-  plot_ly(y = y, z = t(mat_test), x = x, type = "heatmap", colors=palette, showscale=showscale,  text = t(text), hoverinfo='text')%>%
            colorbar(
              thicknessmode="fraction",
              thickness=0.02,
              ypad=25,
              
              tickmode="array", 
              tickvals= tickvals, 
              ticktext= ticktext,
              ticklen=0
            )
          
          
          re_order_seq <- seq(nrow(mat2),1, length=nrow(mat2))
          
          y <- segmented_file[re_order_seq,4]
          
          tickvals <- c(min(as.matrix(round(mat2[re_order_seq,],3))), max_pval)
          range <- tickvals[length(tickvals)]-tickvals[1]
          tickvals[1] <- tickvals[1] + (range*0.1)
          tickvals[length(tickvals)] <- tickvals[length(tickvals)] - (range*0.1)
          
          x <- x
          
          p_pvalues <- plot_ly(y=y, x=x, z = as.matrix(round(mat2[re_order_seq,],3)), type = "heatmap", colors=colfunc, hoverinfo='x+y+z')%>%
            colorbar(
              thicknessmode="fraction",
              thickness=0.02,
              
              title= pval_term,
              titleside = "bottom",
              titlefont=list(size=10, color="black"),
              
              tickmode="array", 
              tickvals=tickvals, 
              ticktext=c(paste("p < ", min(as.matrix(round(mat2[re_order_seq,],3))), sep=""), paste("p > ",max_pval, sep="")),
              ticklen=0
            )%>%
            layout(
              font = list(size=10),
              margin = list(l = 90, r = 90, t = 40, b = 120) 
            )
          
          
          seq_heights <- c(0.2, 0.8)
          p_limma <- subplot(p_groups, p_pvalues, nrows=2, shareX=T, heights=seq_heights, widths=0.9)
          
          output$limma_heatmap <- renderPlotly({
            
            p_limma
            
            ## limma heatmap ##
            #####################################
          })
          
          #download buttons:
          output$dw_buttons_limma_heatmap <- renderUI({
            div(style="text-align:right;",
                downloadButton("TSV_limma_heatmap", "TSV"),
                downloadButton("PNG_limma_heatmap", "PNG")
            )
          })
          
          
          #######################DOWNLOAD HANDLER:####################
          ########TSV_limma_heatmap#####
          output$TSV_limma_heatmap <- downloadHandler(
            filename = paste("p-values_StudentT-test_by_", group_variable, "_groups_in_", segmentation_name, "_", p_val_type, "_", Sys.time(), ".tsv", sep=""),
            content = function(file) {
              write.table(limma_df_2, file, sep="\t", col.names=NA, row.names=T, quote=F, dec=".")
            }
          )
          ########PNG_limma_heatmap#####
          output$PNG_limma_heatmap <- downloadHandler(
            filename = paste("Plot_Student_t-test_by_", group_variable, "_groups_in_", segmentation_name, "_", p_val_type, "_", Sys.time(), ".png", sep=""),
            content = function(file) {
              p_limma$x$layout$font$size <- 70
              p_limma$x$layout$margin$t <- 500
              p_limma$x$layout$margin$b <- 500
              p_limma$x$layout$margin$l <- 500
              p_limma$x$layout$margin$r <- 500
              export(p_limma, file = file, vwidth = 7000, vheight = 5000, expand=1)
              
            }
          )
          
          
          ########################## #
          
        } # else if (length(which(apply(mat2, 2, sd)!=0))>0) (LIMMA)
        
        
        
        #### FISHER'S TEST ####
        
        mat <- data.means
        colnames(mat) <- sample_names
        rownames(mat) <- seg_regions_names
        
        e_list <- list()
        for (i in 1:n_groups){
          name_group <- as.character(real_groups[i])
          samples_g <- which(col_group_var==name_group)
          counts_list <- mclapply(1:nrow(mat), function(x){
            sample_col <- mat[x,samples_g]
            n_gain <- length(which(sample_col >= gain_th)) # gain threshold
            n_loss <- length(which(sample_col <= loss_th)) # loss threshold
            n_wt <- length(which(sample_col < gain_th & sample_col > loss_th))
            ans <- c(n_gain, n_loss, n_wt)
            names(ans) <- c("gain", "loss", "wt")
            ans
          })
          counts <- as.matrix(do.call(rbind, counts_list))
          e_list[[name_group]] <- counts
          
        }
        
        all_counts <- do.call(cbind, e_list)
        rownames(all_counts) <- rownames(mat)
        
        e_vector <- c()
        for (i in 1:length(real_groups)){
          for (j in (i:length(real_groups))[-(i:length(real_groups)==i)]){
            e_vector <- rbind(e_vector, c(real_groups[i], real_groups[j]))
          }
        }
        real_CONTRs <- e_vector
        
        e_vector <- c()
        for (i in 1:length(groups)){
          for (j in (i:length(groups))[-(i:length(groups)==i)]){
            e_vector <- rbind(e_vector, c(groups[i], groups[j]))
          }
        }
        CONTRs <- e_vector
        
        
        
        x <- all_counts
        pval_list_by_real_CONTRs <- list()
        adj.BH_pval_list_by_real_CONTRs <- list()
        for (j in 1:nrow(real_CONTRs)){
          name_vs <- paste(real_CONTRs[j,], collapse="-")
          name.1 <- real_CONTRs[j,1]
          name.2 <- real_CONTRs[j,2]
          cols.1 <- c(1,2,3) + ((which(real_groups==name.1)-1) * 3)
          cols.2 <- c(1,2,3) + ((which(real_groups==name.2)-1) * 3)
          
          pvals <- rep(NA, length(1:nrow(x)))
          for (i in 1:nrow(x)){
            a <- matrix(x[i,cols.1])   
            b <- matrix(x[i,cols.2])
            class <- cbind(b,a)
            rownames(class) <- c("gain", "loss", "wt")
            colnames(class) <- CONTRs[j,]
            
            pvals[i] <- fisher.test(class)$p
            
          }
          p.val.corrected.BH <- p.adjust(pvals, method="BH")
          pval_list_by_real_CONTRs[[name_vs]] <- pvals
          adj.BH_pval_list_by_real_CONTRs[[name_vs]] <- p.val.corrected.BH
          
          
          
        }
        
        
        if (p_val_type == "p-value") {
          pval_term <- "p-val"
          
          ans_fisher <- as.matrix(do.call(cbind, pval_list_by_real_CONTRs))
          rownames(ans_fisher) <- segmented_file[,4]
        } else if (p_val_type == "Adj.p-value") {
          pval_term <- "Adj.p-val"
          
          ans_fisher <- as.matrix(do.call(cbind, adj.BH_pval_list_by_real_CONTRs))
          rownames(ans_fisher) <- segmented_file[,4]
        }
        
        annot_fisher <- as.matrix(real_CONTRs)
        
        e_vector <- rep(NA, length=1)
        for (i in 1:length(groups)){
          for (j in (i:length(groups))[-(i:length(groups)==i)]){
            contr <- paste(groups[j],groups[i],sep="-")
            e_vector <- c(e_vector, contr)
          }
        }
        collapse_CONTRs <- e_vector[-1]
        
        rownames(annot_fisher) <- collapse_CONTRs
        colnames(annot_fisher) <- c("Group.1","Group.2")
        
        
        
        
        ### output ####
        
        output$fisher_heatmap_MAIN <- renderUI({
          HTML(paste("<br><h4>More differentiated regions between groups defined by <i><b>",
                     group_variable, "</b></i>", sep=""))
          
        })
        
        
        n.colors <- length(seq(by=0.01, from=-2, to=2)) # centrat a 0
        #	n.colors <- round((MAX-MIN)+1)+100 # centrat a la mediana de la sèrie entre MAX i MIN
        colfunc <- colorRampPalette(c("red", "yellow", "white"))(n.colors)
        
        mat2 <- as.matrix(ans_fisher)
        #rownames(mat2) <- paste(segmented_file[,1], segmented_file[,4])
        for (i in 1:nrow(mat2)) {
          for (j in 1:ncol(mat2)) {
            #if (mat2[i,j]<(0.05)) mat2[i,j]<- (0.05)
            if (is.na(mat2[i,j])) mat2[i,j]<- 1
            
            if (mat2[i,j]>(max_pval)) mat2[i,j]<- (max_pval)
            if (mat2[i,j]<(min_pval)) mat2[i,j]<- (min_pval)
            
          }
        }
        
        if (length(which(apply(mat2, 2, sd)!=0))==0){
          #output$fihser_heatmap_MAIN <- renderUI({HTML("")})
          #output$text_limma <- renderUI(div(style="color:red","Try to expand p-value range!!"))
          output$not_differences_2_fisher <- output$not_differences_fisher <- renderUI(HTML(paste("<p style='color:red'>Sorry, but there are no descriptive regions between '", group_variable, "' groups<br>Try to expand p-value range!!</p>")))
          output$fisher_heatmap <- renderPlotly({
            
            plotly_empty()
            
          })
          output$dw_buttons_fisher_heatmap <- renderUI({
            div(style="text-align:right;",
                actionButton("empty_1", "TSV"),
                actionButton("empty_2", "PNG")
            )
          })
          # output$dw_buttons_reg_boxplot <- renderUI({
          #   div(style="text-align:right;",
          #       actionButton("empty_3", "TSV"),
          #       actionButton("empty_4", "PNG")
          #   )
          # })
          
          
          
        } else if (length(which(apply(mat2, 2, sd)!=0))>0) {
          output$not_differences_2_fisher <- output$not_differences_fisher <- renderUI(div(""))
          
          ####the plot####
          
          
          m <- as.matrix(t(annot_fisher))
          m1 <- mapply(m, FUN=as.factor)
          data.groups <- matrix(data=as.numeric(m1), ncol=ncol(m), nrow=nrow(m))
          colnames(data.groups) <- rownames(annot_fisher)
          rownames(data.groups) <- colnames(annot_fisher)
          
          y <- rownames(data.groups)
          palette <- "Spectral"
          palette <- as.character(variables_info_mat[which(variables_info_mat[,"name_var"]==group_variable),"color_palette"])
          showscale <- TRUE
          
          x <- paste("'", names(ans_fisher), "'", sep="")
          
          text <- t(m)
          
          mat_test <- t(data.groups)
          tickvals <- unique(as.numeric(m1))
          range <- tickvals[length(tickvals)]-tickvals[1]
          tickvals[1] <- tickvals[1] + (range*0.1)
          tickvals[length(tickvals)] <- tickvals[length(tickvals)] - (range*0.1)
          ticktext <- levels(m1)
          
          
          p_groups <-  plot_ly(y = y, z = t(mat_test), x = x, type = "heatmap", colors=palette, showscale=showscale,  text = t(text), hoverinfo='text')%>%
            colorbar(
              thicknessmode="fraction",
              thickness=0.02,
              ypad=25,
              
              tickmode="array", 
              tickvals= tickvals, 
              ticktext= ticktext,
              ticklen=0
            )
          
          
          re_order_seq <- seq(nrow(mat2),1, length=nrow(mat2))
          
          y <- segmented_file[re_order_seq,4]
          x <- x
          
          tickvals <- c(min(as.matrix(round(mat2[re_order_seq,],3))), max_pval)
          range <- tickvals[length(tickvals)]-tickvals[1]
          tickvals[1] <- tickvals[1] + (range*0.1)
          tickvals[length(tickvals)] <- tickvals[length(tickvals)] - (range*0.1)
          
          
          p_pvalues <- plot_ly(y=y, x=x, z = as.matrix(round(mat2[re_order_seq,],3)), type = "heatmap", colors=colfunc, hoverinfo='x+y+z')%>%
            colorbar(
              thicknessmode="fraction",
              thickness=0.02,
              
              title= pval_term,
              titleside = "bottom",
              titlefont=list(size=10, color="black"),
              
              tickmode="array", 
              tickvals=tickvals, 
              ticktext=c(paste("p < ", min(as.matrix(round(mat2[re_order_seq,],3))), sep=""), paste("p > ",max_pval, sep="")),
              ticklen=0
            )%>%
            layout(
              font = list(size=10),
              margin = list(l = 90, r = 90, t = 40, b = 120) 
            )
          
          
          seq_heights <- c(0.2, 0.8)
          p_fisher <- subplot(p_groups, p_pvalues, nrows=2, shareX=T, heights=seq_heights, widths=0.9)
          
          output$fisher_heatmap <- renderPlotly({
            
            p_fisher
            
          })
          
          #download buttons:
          output$dw_buttons_fisher_heatmap <- renderUI({
            div(style="text-align:right;",
                downloadButton("TSV_fisher_heatmap", "TSV"),
                downloadButton("PNG_fisher_heatmap", "PNG")
            )
          })
          
          
          #######################DOWNLOAD HANDLER:####################
          ########TSV_fisher_heatmap#####
          output$TSV_fisher_heatmap <- downloadHandler(
            filename = paste("p-values_Fisher_t-test_by_", group_variable, "_groups_in_", segmentation_name, "_", p_val_type, "_", Sys.time(), ".tsv", sep=""),
            content = function(file) {
              write.table(ans_fisher, file, sep="\t", col.names=NA, row.names=T, quote=F, dec=".")
            }
          )
          ########PNG_fisher_heatmap#####
          output$PNG_fisher_heatmap <- downloadHandler(
            filename = paste("Plot_Fisher_t-test_by_", group_variable, "_groups_in_", segmentation_name, "_", p_val_type, "_", Sys.time(), ".png", sep=""),
            content = function(file) {
              p_fisher$x$layout$font$size <- 70
              p_fisher$x$layout$margin$t <- 500
              p_fisher$x$layout$margin$b <- 500
              p_fisher$x$layout$margin$l <- 500
              p_fisher$x$layout$margin$r <- 500
              export(p_fisher, file = file, vwidth = 7000, vheight = 5000, expand=1)
              
            })
          
          
          ########################## #
          
        } # else if (length(which(apply(mat2, 2, sd)!=0))>0)  (FISHER)
        
        ############ OUTPUT Boxplot by region############
        p_boxplot <- reactive({
          region <- as.character(input$reg_selected)
          region
          
          data_2 <- as.data.frame(t(rbind(data.means,as.character(mat_variables[,group_variable]))))
          colnames(data_2) <- c(seg_regions_names,group_variable)
          
          data_3 <- data_2[which(data_2[,group_variable]%in%real_groups),]
          for (i in 1:(ncol(data_3)-1)){
            data_3[,i] <- as.numeric(as.character(data_3[,i]))
          }
          
          if (length(real_groups)<=6) {
            
            e_list <- list()
            for (i in 1:length(real_groups)){
              for (j in (i:length(real_groups))[-(i:length(real_groups)==i)]){
                contr <- c(as.character(real_groups[i]),as.character(real_groups[j]))
                e_list <- c(e_list, list(contr))
              }
            }
            
            p_vals_here <- round(as.numeric(as.character(limma_df[which(rownames(limma_df)==as.character(region)),])),3)
            
            for (i in 1:length(p_vals_here)){
              val <- as.numeric(as.character(p_vals_here[i]))
              if (val<=0.0001){
                val2 <- paste("****", val, sep=" ")
              } else if (val <= 0.001){
                val2 <- paste("***", val, sep=" ")
              } else if (val <= 0.01){
                val2 <- paste("**", val, sep=" ")
              } else if (val <= 0.05){
                val2 <- paste("*", val, sep=" ")
              } else if (val > 0.05){
                val2 <- paste("(ns)", val, sep=" ")
              } else {
                val2 <- val
              }
              p_vals_here[i] <- val2
            }
            
            p <- ggplot(data_3, aes(x=reorder(get(group_variable)), y=get(reg_to_plot()))) + geom_boxplot()
            p_boxplot <- p + xlab(group_variable) + ylab(as.character(region)) + 
              geom_signif(comparisons = e_list, map_signif_level=FALSE, step_increase=0.08 , annotation=p_vals_here) +
              theme(text = element_text(size=14), axis.text.x = element_text(size=12), plot.margin = unit(c(1,1,1,1), "cm"))
            
          } else if (length(real_groups)>6) {
            p <- ggplot(data_3, aes(x=reorder(get(group_variable)), y=get(reg_to_plot()))) + geom_boxplot()
            p_boxplot <- p + xlab(group_variable) + ylab(as.character(region)) + 
              theme(text = element_text(size=14), axis.text.x = element_text(size=12), plot.margin = unit(c(1,1,1,1), "cm"))
          }
          
          return(p_boxplot)
          
        })
        
        output$reg_boxplot_MAIN_fisher <- output$reg_boxplot_MAIN <- renderUI({
          HTML(paste("<br><h4>Boxplot for region <u><b>",
                     as.character(reg_to_plot()),
                     "</b></u> (<i style='color:orange'>by ", segmentation_name,"</i>)</h4>", sep=""))
        })
        
        
        
        output$reg_boxplot <- renderPlot({
          p_boxplot()
        })
        
        #download buttons boxplots:
        output$dw_buttons_reg_boxplot_fisher <- output$dw_buttons_reg_boxplot <- renderUI({
          div(style="text-align:right;",
              downloadButton("TSV_reg_boxplot", "TSV"),
              downloadButton("PNG_reg_boxplot", "PNG")
          )
        })
        
        ############ OUTPUT Stacked by region############
        p_stacked <- reactive({
          region <- as.character(input$reg_selected)
          region
          
          mt.0 <- all_counts[which(rownames(all_counts)==region),]
          fun_here <- function(i) {
            name_group <- real_groups[i]
            j <- 3*(i-1)+1
            t <- mt.0[j:(j+2)]
            sum_t <- sum(t)
            ans <- round((t/sum_t)*100,1)
          }
          mt.1 <- t(as.matrix(sapply(1:length(real_groups), fun_here)))
          mt <- as.data.frame(cbind(real_groups, mt.1))
          rownames(mt) <- real_groups
          
          mt$real_groups <- factor(mt$real_groups, levels = c(as.character(mt$real_groups)))
          
          
          colors <- c("red", "blue", "gainsboro")
          #par(xpd=T, mar=c(1,3,1,5))
          
          x <- as.character(rownames(mt))
          
          p_stacked <- plot_ly(mt, x=~real_groups, y=mt[,"gain"], name="gain",type = "bar", marker=list(color=colors[1]), showlegend=F)%>%
            add_trace(y=mt[,"loss"], name="loss", marker=list(color=colors[2]))%>%
            add_trace(y=mt[,"wt"], name="wt", marker=list(color=colors[3]))%>%
            
            layout(
              yaxis =list(title= "(%) CNA events", type="linear", range=seq(1:sum(apply(mt.1,1,sum)))),
              xaxis=list(title=group_variable),
              barmode="stack",
              
              font = list(size=10),
              margin = list(l = 90, r = 90, t = 80, b = 120)
            )
          return(p_stacked)
        })
        
        output$reg_stacked_MAIN <- renderUI({
          HTML(paste("<br><h4>% of CNA events in region <u><b>",
                     as.character(reg_to_plot()),
                     "</b></u> (<i style='color:orange'>by ", segmentation_name,"</i>)</h4>", sep=""))
        })
        
        
        
        output$reg_stacked <- renderPlotly({
          p_stacked()
        })
        
        #download button stacked:
        output$dw_buttons_reg_stacked <- renderUI({
          div(style="text-align:right;",
              downloadButton("TSV_reg_stacked", "TSV"),
              downloadButton("PNG_reg_stacked", "PNG")
          )
        })
        
        
        #### Data genes in region
        
        t_genes <- reactive({
          term <- as.character(input$reg_selected)
          term
          
          row_region <- segmented_file[which(segmented_file[,4]==term),]
          
          fun_name <- "refGene.selection" # function name
          my_fun <- paste(fun_name, ".R", sep="")
          source_fun <- paste(dir_funs, "/", my_fun, sep="")
          source(source_fun)
          # sourcing fun into R
          
          refGene_mat <- refGene.selection(genome_build)
          
          CHROM <- as.character(row_region[,1])
          START <- as.numeric(as.character(row_region[,2]))
          END <- as.numeric(as.character(row_region[,3]))
          
          rows_in_refGene <- which(refGene_mat[,"chrom"]==CHROM & as.numeric(as.character(refGene_mat[,"start"]))>=START & as.numeric(as.character(refGene_mat[,"end"]))<=END)
          
          
          
          genes_in_reg <- refGene_mat[rows_in_refGene,c("symbol_name", "chrom", "strand", "start", "end")]
          uni_g <- unique(as.character(genes_in_reg[,1]))
          each_gene <- as.vector(do.call(cbind, mclapply(uni_g, function(x,y){which(y==x)[1]}, y=as.character(genes_in_reg[,1]))))
          
          genes_in_reg_2 <- genes_in_reg[each_gene,]
          genes_in_reg_2
        })
        
        output$t_genes_reg_MAIN_fisher <- output$t_genes_reg_MAIN <- renderUI({
          HTML(paste("<br><h4>Genes in region <u><b>", as.character(reg_to_plot()), sep=""))
        })
        
        output$t_genes_reg_fisher <- output$t_genes_reg <- renderDataTable({
          return(t_genes())
        })
        
        #download buttons boxplots:
        output$dw_tab_genes_in_reg_fisher <- output$dw_tab_genes_in_reg <- renderUI({
          div(style="text-align:right;",
              downloadButton("TSV_tab_genes_in_reg", "TSV")
          )
        })
        
        
        #######Download reg_boxplot######
        ########PNG_reg_boxplot#####
        output$PNG_reg_boxplot <- downloadHandler(
          filename = reactive({
            term <- as.character(input$reg_selected)
            paste("Region_boxplot_in_", group_variable, "_groups_", as.character(term), "_", p_val_type, "_", segmentation_name, Sys.time(), ".png", sep="")
          }),
          content = function(file) {
            # p_boxplot[["theme"]][["text"]][["size"]] <- 70
            # p_boxplot[["theme"]][["axis.text.x"]][["size"]] <- 65
            # p_boxplot[["theme"]][["plot.margin"]] <- c(2,2,2,2)
            
            p_boxplot <- reactive({
              region <- as.character(input$reg_selected)
              region
              
              data_2 <- as.data.frame(t(rbind(data.means,as.character(mat_variables[,group_variable]))))
              colnames(data_2) <- c(seg_regions_names,group_variable)
              
              data_3 <- data_2[which(data_2[,group_variable]%in%real_groups),]
              for (i in 1:(ncol(data_3)-1)){
                data_3[,i] <- as.numeric(as.character(data_3[,i]))
              }
              
              if (length(real_groups)<=6) {
                
                e_list <- list()
                for (i in 1:length(real_groups)){
                  for (j in (i:length(real_groups))[-(i:length(real_groups)==i)]){
                    contr <- c(as.character(real_groups[i]),as.character(real_groups[j]))
                    e_list <- c(e_list, list(contr))
                  }
                }
                
                p_vals_here <- round(as.numeric(as.character(limma_df[which(rownames(limma_df)==as.character(region)),])),3)
                
                for (i in 1:length(p_vals_here)){
                  val <- as.numeric(as.character(p_vals_here[i]))
                  if (val<=0.0001){
                    val2 <- paste("****", val, sep=" ")
                  } else if (val <= 0.001){
                    val2 <- paste("***", val, sep=" ")
                  } else if (val <= 0.01){
                    val2 <- paste("**", val, sep=" ")
                  } else if (val <= 0.05){
                    val2 <- paste("*", val, sep=" ")
                  } else if (val > 0.05){
                    val2 <- paste("(ns)", val, sep=" ")
                  } else {
                    val2 <- val
                  }
                  p_vals_here[i] <- val2
                }
                
                p <- ggplot(data_3, aes(x=get(group_variable), y=get(reg_to_plot()))) + geom_boxplot()
                p_boxplot <- p + xlab(group_variable) + ylab(as.character(region)) + 
                  geom_signif(comparisons = e_list, map_signif_level=FALSE, step_increase=0.08 , annotation=p_vals_here, size=1, textsize = 8) +
                  theme(text = element_text(size=30), axis.text.x = element_text(size=30), plot.margin = unit(c(1,1,1,1), "cm"))
                
              } else if (length(real_groups)>6) {
                p <- ggplot(data_3, aes(x=get(group_variable), y=get(reg_to_plot()))) + geom_boxplot()
                p_boxplot <- p + xlab(group_variable) + ylab(as.character(region)) + 
                  theme(text = element_text(size=30), axis.text.x = element_text(size=30), plot.margin = unit(c(1,1,1,1), "cm"))
              }
              
              return(p_boxplot)
              
            })
            
            png(file, width = 7000, height = 5000, res = 400)
            plot(p_boxplot())
            dev.off()
            #plotPNG(func = plot(p_boxplot()), filename = file, width = 7000, height = 5000, res=500)
            
          }
        )
        ########TSV_reg_boxplot#####
        output$TSV_reg_boxplot <- downloadHandler(
          filename = reactive({
            term <- as.character(input$reg_selected)
            paste("Region_data_in_", group_variable, "_groups_", as.character(term), "_", p_val_type, "_", segmentation_name, Sys.time(), ".tsv", sep="")
          }),
          content = function(file) {
            
            boxplot_data <- reactive({
              region <- as.character(input$reg_selected)
              region
              
              data_2 <- as.data.frame(t(rbind(data.means,mat_variables[,group_variable])))
              colnames(data_2) <- c(seg_regions_names,group_variable)
              
              data_3 <- data_2[which(data_2[,group_variable]%in%real_groups),]
              for (i in 1:(ncol(data_3)-1)){
                data_3[,i] <- as.numeric(as.character(data_3[,i]))
              }
              data_4 <- as.data.frame(cbind(colnames(data.means), data_3[,c(region, group_variable)]))
              colnames(data_4) <- c("ID", region, group_variable)
              return(data_4)
              
            })
            
            region <- as.character(input$reg_selected)
            region
            
            data_2 <- as.data.frame(t(rbind(data.means,mat_variables[,group_variable])))
            colnames(data_2) <- c(seg_regions_names,group_variable)
            
            data_3 <- data_2[which(data_2[,group_variable]%in%real_groups),]
            for (i in 1:(ncol(data_3)-1)){
              data_3[,i] <- as.numeric(as.character(data_3[,i]))
            }
            
            write.table(as.data.frame(boxplot_data()), file, sep="\t", col.names=T, row.names=F, quote=F, dec=".")
          }
        )
        
        
        ########################## #
        ########PNG_reg_stacked#####
        output$PNG_reg_stacked <- downloadHandler(
          filename = reactive({
            term <- as.character(input$reg_selected)
            paste("Region_fisher_", group_variable, "_groups_", as.character(term), "_", p_val_type, "_", segmentation_name, Sys.time(), ".png", sep="")
          }),
          content = function(file) {
            # p_boxplot[["theme"]][["text"]][["size"]] <- 70
            # p_boxplot[["theme"]][["axis.text.x"]][["size"]] <- 65
            # p_boxplot[["theme"]][["plot.margin"]] <- c(2,2,2,2)
            
            
            #export(p_boxplot, file = file, vwidth = 7000, vheight = 5000, expand=1)
            p_stacked <- p_stacked()
            p_stacked$x$layout$font$size <- 200
            p_stacked$x$layout$margin$t <- 700
            p_stacked$x$layout$margin$b <- 1100
            p_stacked$x$layout$margin$l <- 1000
            p_stacked$x$layout$margin$r <- 700
            export(p_stacked,file = file, vwidth = 1200, vheight = 700, expand=1)
            
            
          }
        )
        ########TSV_reg_stacked#####
        output$TSV_reg_stacked <- downloadHandler(
          filename = reactive({
            term <- as.character(input$reg_selected)
            paste("Region_fisher_data_", group_variable, "_groups_", as.character(term), "_", p_val_type, "_", segmentation_name, Sys.time(), ".tsv", sep="")
          }),
          content = function(file) {
            # p_boxplot[["theme"]][["text"]][["size"]] <- 70
            # p_boxplot[["theme"]][["axis.text.x"]][["size"]] <- 65
            # p_boxplot[["theme"]][["plot.margin"]] <- c(2,2,2,2)
            
            reg_counts <- reactive({
              region <- as.character(input$reg_selected)
              region
              
              mt.0 <- all_counts[which(rownames(all_counts)==region),]
              fun_here <- function(i) {
                name_group <- real_groups[i]
                j <- 3*(i-1)+1
                t <- mt.0[j:(j+2)]
                t
                # sum_t <- sum(t)
                # ans <- round((t/sum_t)*100,1)
              }
              mt.1 <- t(as.matrix(sapply(1:length(real_groups), fun_here)))
              mt <- as.data.frame(cbind(real_groups, mt.1))
              rownames(mt) <- real_groups
              
              
              return(mt)
            })
            
            write.table(reg_counts(), file, sep="\t", col.names=T, row.names=F, quote=F, dec=".")
            
          }
        )
        
        ########################## #
        ########TSV_tab_genes_in_reg#####
        reg_name <- reactive({
          term <- as.character(input$reg_selected)
          paste("Genes_in_reg_", as.character(term), "_", segmentation_name, Sys.time(), ".tsv", sep="")
        })
        
        output$TSV_tab_genes_in_reg_fisher <- output$TSV_tab_genes_in_reg <- downloadHandler(
          filename = reactive({
            term <- as.character(input$reg_selected)
            paste("Genes_in_reg_", as.character(term), "_", segmentation_name, Sys.time(), ".tsv", sep="")
          }),
          content = function(file) {
            write.table(as.data.frame(t_genes()), file, sep="\t", col.names=T, row.names=F, quote=F, dec=".")
          }
        )
        
      }#else if (n_groups > 1)
      
      
      
      
      
      
    })#observeEvent (button_run_limma)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    all_profile_part <<- profile_part
    
    output$run_parts_in_create_model_profile <- renderUI({
      div(
        radioButtons("profile_param_for_model", label=HTML("Activate <a style=background-color:white>&nbsp<i class='fa fa-map' style='color:orange'></i>&nbsp <i style=color:gray>REGION PROFILE &nbsp</i></a> parameters:"), 
                     choices=c("Yes", "No"), 
                     selected="Yes",
                     inline = TRUE),
        if (q_use_cin_ranges=="Yes") {
          specific_ranges_side <- specific_ranges
          if (specific_ranges_side=="Broad (chromosomal and arm-level)") { specific_ranges_side <- "Broad" }
          HTML("<h5 style=text-align:left;color:black>Using <i><b style=color:green>", specific_ranges_side, " CNAs</b></i> from &nbsp<i class='fa fa-star' style='color:orange'></i>&nbsp <i style=color:gray>RE-SEG & SCORE &nbsp</i></h5>
           <br>")
        } else if (q_use_cin_ranges=="No") {
          HTML("")
        }
      )
    })
    output$run_parts_in_create_model_data_loaded <- renderUI({
      div(
        br(),
        div(style="text-align:left", title="click 'Load' to load new parameters", actionBttn("button_param_create_model", "Load", size="sm", color="primary", style="simple"))
      )
    })
    
  })#observeEvent (button_run)
  
  
  ############################
  ############################
  ########### CREATE YOUR MODEL param#################
  
  observeEvent(input$button_param_create_model, { #when 'Load' is clicked!
    
    output$global_df <- renderTable(GLOBAL_DF)
    
    
    use_data_loaded <<- "Yes"
    
    if (cin_part_run == "YES" & cn_profile_part_run == "YES") {
      
      use_cin_param <- as.character(input$cin_param_for_model)
      use_profile_param <- as.character(input$profile_param_for_model)
      
      
    } else if (cin_part_run == "YES" & cn_profile_part_run == "NO") {
      
      use_cin_param <- as.character(input$cin_param_for_model)
      use_profile_param <- "No"
      
      
    } else if (cin_part_run == "NO" & cn_profile_part_run == "YES") {
      
      use_cin_param <- "No"
      use_profile_param <- as.character(input$profile_param_for_model)
      
      
    } else if (cin_part_run == "NO" & cn_profile_part_run == "NO") {
      
      use_cin_param <- "No"
      use_profile_param <- "No"
      
      
    }
    
    
    vars_option <<- "READ"
    
    if (use_cin_param == "Yes" & use_profile_param == "Yes" & use_data_loaded == "Yes") {
      vars_option <<- "CIN_PROFILE_READ"
      
      vars_cin <- colnames(all_cin_part$scores)
      vars_profile <- colnames(all_profile_part$mat_variables)
      repeated_vars <- which(vars_profile%in%vars_cin)
      if (length(repeated_vars)>0){
        vars_profile <- vars_profile[-repeated_vars]
      }
      vars_profile <- vars_profile[-which(vars_profile=="ID")]
      
      sel_list_1 <- list("Annotation variables" = vars_profile, "CIN SCOREs" = vars_cin)
      
      sel_list_2 <- list("Regions from profile" = c("All regions", "Specific region/s"), "CIN SCOREs" = vars_cin, "Annotation variables" = vars_profile)
      
      output$create_model_param <- renderUI({
        div(
          HTML("<h4>Choose your variables:</h4>"),
          
          fluidRow(
            column(6,
                   #variable to define sample groups in model
                   selectizeInput("var_def_groups", "Variable to DEFINE groups", choices= sel_list_1, multiple=FALSE, options = list(maxOptions=3000))
            ),
            column(6,
                   #variable (or variables...) to classify sample groups defined by 'var_def_groups'
                   selectizeInput("var_to_classify", "Variable/s to CLASSIFY groups", choices= sel_list_2, multiple=TRUE, selected= sel_list_2[[1]][1], options = list(maxOptions=3000))
            )
          )
        )
      })
      
      subs_variable_group <- eventReactive(input$var_def_groups, {
        groups_by_var <- as.character(input$var_def_groups)
        if (length(grep(groups_by_var, vars_profile))>0  ) {
          where <- "profile"
          class_var <- as.character(all_profile_part$variables_info_mat[which(all_profile_part$variables_info_mat[,1]==groups_by_var),2])
          if (class_var=="categoric") {
            legend_var <- sort(unique(as.character(all_profile_part$mat_variables[,colnames(all_profile_part$mat_variables)==groups_by_var])))
            div(
              HTML(paste("<i>",length(legend_var), " groups defined:</i>", sep="")),
              checkboxGroupInput("sub_group_var", "", choices = legend_var, selected = legend_var)
            )
          } else if (class_var=="numeric"){
            values <- round(as.numeric(as.character(all_profile_part$mat_variables[,colnames(all_profile_part$mat_variables)==groups_by_var])),3)
            div(
              HTML(paste("<i>",length(unique(values)), " groups defined:</i>", sep="")),
              sliderInput("slider_group_var_model", "", min=min(values, na.rm = T), max=max(values, na.rm = T), value = c(min(values, na.rm = T),max(values, na.rm = T))),
              sliderInput("slider_group_var_model_2", "", min=min(values, na.rm = T), max=max(values, na.rm = T), value = c(min(values, na.rm = T),max(values, na.rm = T))) 
            )
          }
          
        } else if (length(grep(groups_by_var, vars_cin))>0  ) {
          where <- "cin"
          values <- round(as.numeric(as.character(all_cin_part$scores[,colnames(all_cin_part$scores)==groups_by_var])),3)
          div(
            HTML(paste("<i>",length(unique(values)), " groups defined:</i>", sep="")),
            sliderInput("slider_group_var_model", "", min=min(values, na.rm = T), max=max(values, na.rm = T), value = c(min(values, na.rm = T),max(values, na.rm = T))),
            sliderInput("slider_group_var_model_2", "", min=min(values, na.rm = T), max=max(values, na.rm = T), value = c(min(values, na.rm = T),max(values, na.rm = T))) 
          )
          
        }
        
        
      })
      
      
      spec_out <- eventReactive(input$var_to_classify, {
        to_classify_var <- input$var_to_classify
        
        if (length( grep("Specific region/s", to_classify_var) ) ==1 ) {
          
          selectizeInput("spec_regions", "Select specific regions:", choices=as.character(all_profile_part[["segmented_file"]][,4]), selected=as.character(all_profile_part[["segmented_file"]][1:3,4]), multiple=TRUE, options = list(maxOptions=3000))
        } else {
          
          HTML("")
          
        }
        
      })
      
      
      output$create_model_param_2 <- renderUI({
        fluidRow(
          column(6,
                 subs_variable_group(),
                 div(style="text-align:left",actionBttn("button_create_model", "Create new model", style="fill", size="md", color="warning")),
                 busyIndicator(text="Running", wait=200)
          ),
          column(6,
                 spec_out()
          )
        )
        
      })
      
      
    } else if (use_cin_param == "Yes" & use_profile_param == "No" & use_data_loaded == "Yes") {
      vars_option <<- "CIN_READ"
      
      vars_cin <- colnames(all_cin_part$scores)
      annot_vars <-  colnames(all_read_part$mat_variables.0)
      repeated_vars <- which(annot_vars%in%vars_cin)
      if (length(repeated_vars)>0){
        annot_vars <- annot_vars[-repeated_vars]
      }
      annot_vars <- annot_vars[-which(annot_vars=="ID")]
      
      sel_list_1 <- list("Annotation variables" = annot_vars, "CIN SCOREs" = vars_cin)
      
      sel_list_2 <- list("CIN SCOREs" = vars_cin, "Annotation variables" = annot_vars)
      
      
      output$create_model_param <- renderUI({
        div(
          HTML("<h4>Choose your variables:</h4>"),
          
          fluidRow(
            column(6,
                   #variable to define sample groups in model
                   selectizeInput("var_def_groups", "Variable to DEFINE groups", choices= sel_list_1, multiple=FALSE, options = list(maxOptions=3000))
            ),
            column(6,
                   #variable (or variables...) to classify sample groups defined by 'var_def_groups'
                   selectizeInput("var_to_classify", "Variable/s to CLASSIFY groups", choices= sel_list_2, multiple=TRUE, selected= sel_list_2[[1]][1], options = list(maxOptions=3000))
            )
          )
        )
      })
      
      subs_variable_group <- eventReactive(input$var_def_groups, {
        groups_by_var <- as.character(input$var_def_groups)
        if (length(grep(groups_by_var, annot_vars))>0  ) {
          where <- "read"
          class_var <- as.character(all_read_part$variables_info_mat[which(all_read_part$variables_info_mat[,1]==groups_by_var),2])
          if (class_var=="categoric") {
            legend_var <- sort(unique(as.character(all_read_part$mat_variables.0[,colnames(all_read_part$mat_variables.0)==groups_by_var])))
            div(
              HTML(paste("<i>",length(legend_var), " groups defined:</i>", sep="")),
              checkboxGroupInput("sub_group_var", "", choices = legend_var, selected = legend_var)
            )
          } else if (class_var=="numeric"){
            values <- round(as.numeric(as.character(all_read_part$mat_variables.0[,colnames(all_read_part$mat_variables.0)==groups_by_var])),3)
            div(
              HTML(paste("<i>",length(unique(values)), " groups defined:</i>", sep="")),
              sliderInput("slider_group_var_model", "", min=min(values, na.rm = T), max=max(values, na.rm = T), value = c(min(values, na.rm = T),max(values, na.rm = T))),
              sliderInput("slider_group_var_model_2", "", min=min(values, na.rm = T), max=max(values, na.rm = T), value = c(min(values, na.rm = T),max(values, na.rm = T))) 
            )
          }
          
        } else if (length(grep(groups_by_var, vars_cin))>0  ) {
          where <- "cin"
          values <- round(as.numeric(as.character(all_cin_part$scores[,colnames(all_cin_part$scores)==groups_by_var])),3)
          div(
            HTML(paste("<i>",length(unique(values)), " groups defined:</i>", sep="")),
            sliderInput("slider_group_var_model", "", min=min(values, na.rm = T), max=max(values, na.rm = T), value = c(min(values, na.rm = T),max(values, na.rm = T))),
            sliderInput("slider_group_var_model_2", "", min=min(values, na.rm = T), max=max(values, na.rm = T), value = c(min(values, na.rm = T),max(values, na.rm = T))) 
          )
          
        }
        
        
      })
      
      output$create_model_param_2 <- renderUI({
        fluidRow(
          column(6,
                 subs_variable_group(),
                 div(style="text-align:left",actionBttn("button_create_model", "Create new model", style="fill", size="md", color="warning")),
                 busyIndicator(text="Running", wait=200)
          ),
          column(6
                 
          )
        )
        
      })
      
      
      
    } else if (use_cin_param == "No" & use_profile_param == "Yes" & use_data_loaded == "Yes"){
      vars_option <<- "PROFILE_READ"
      
      vars_profile <- colnames(all_profile_part$mat_variables.0)
      annot_vars <-  colnames(all_read_part$mat_variables.0)
      repeated_vars <- which(annot_vars%in%vars_profile)
      if (length(repeated_vars)>0){
        annot_vars <- annot_vars[-repeated_vars]
      }
      annot_vars <- annot_vars[-which(annot_vars=="ID")]
      
      sel_list_1 <- list("Annotation variables" = c(vars_profile, annot_vars))
      
      sel_list_2 <- list("Regions from profile" = c("All regions", "Specific region/s"), "Annotation variables" = c(vars_profile, annot_vars))
      
      
      output$create_model_param <- renderUI({
        div(
          HTML("<h4>Choose your variables:</h4>"),
          
          fluidRow(
            column(6,
                   #variable to define sample groups in model
                   selectizeInput("var_def_groups", "Variable to DEFINE groups", choices= sel_list_1, multiple=FALSE, options = list(maxOptions=3000))
            ),
            column(6,
                   #variable (or variables...) to classify sample groups defined by 'var_def_groups'
                   selectizeInput("var_to_classify", "Variable/s to CLASSIFY groups", choices= sel_list_2, multiple=TRUE, selected= sel_list_2[[1]][1], options = list(maxOptions=3000))
            )
          )
        )
      })
      
      subs_variable_group <- eventReactive(input$var_def_groups, {
        groups_by_var <- as.character(input$var_def_groups)
        if (length(grep(groups_by_var, annot_vars))>0  ) {
          where <- "read"
          class_var <- as.character(all_read_part$variables_info_mat[which(all_read_part$variables_info_mat[,1]==groups_by_var),2])
          if (class_var=="categoric") {
            legend_var <- sort(unique(as.character(all_read_part$mat_variables.0[,colnames(all_read_part$mat_variables.0)==groups_by_var])))
            div(
              HTML(paste("<i>",length(legend_var), " groups defined:</i>", sep="")),
              checkboxGroupInput("sub_group_var", "", choices = legend_var, selected = legend_var)
            )
          } else if (class_var=="numeric"){
            values <- round(as.numeric(as.character(all_read_part$mat_variables.0[,colnames(all_read_part$mat_variables.0)==groups_by_var])),3)
            div(
              HTML(paste("<i>",length(unique(values)), " groups defined:</i>", sep="")),
              br(),
              sliderInput("slider_group_var_model", "", min=min(values, na.rm = T), max=max(values, na.rm = T), value = c(min(values, na.rm = T),max(values, na.rm = T))),
              sliderInput("slider_group_var_model_2", "", min=min(values, na.rm = T), max=max(values, na.rm = T), value = c(min(values, na.rm = T),max(values, na.rm = T))) 
            )
          }
          
        } else if (length(grep(groups_by_var, vars_profile))>0  ) {
          where <- "profile"
          values <- round(as.numeric(as.character(all_profile_part$mat_variables.0[,colnames(all_profile_part$mat_variables.0)==groups_by_var])),3)
          div(
            HTML(paste("<i>",length(unique(values)), " groups defined:</i>", sep="")),
            sliderInput("slider_group_var_model", "", min=min(values, na.rm = T), max=max(values, na.rm = T), value = c(min(values, na.rm = T),max(values, na.rm = T))),
            sliderInput("slider_group_var_model_2", "", min=min(values, na.rm = T), max=max(values, na.rm = T), value = c(min(values, na.rm = T),max(values, na.rm = T))) 
          )
          
        }
        
      })
      
      spec_out <- eventReactive(input$var_to_classify, {
        to_classify_var <- input$var_to_classify
        
        if (length( grep("Specific region/s", to_classify_var) ) ==1 ) {
          
          selectizeInput("spec_regions", "Select specific regions:", choices=as.character(all_profile_part[["segmented_file"]][,4]), selected=as.character(all_profile_part[["segmented_file"]][1:3,4]), multiple=TRUE, options = list(maxOptions=3000))
        } else {
          
          HTML("")
          
        }
        
      })
      
      
      output$create_model_param_2 <- renderUI({
        fluidRow(
          column(6,
                 subs_variable_group(),
                 div(style="text-align:left",actionBttn("button_create_model", "Create new model", style="fill", size="md", color="warning")),
                 busyIndicator(text="Running", wait=200)
          ),
          column(6,
                 spec_out()
          )
        )
        
      })
      
      
      
      
    } else if (use_cin_param == "No" & use_profile_param == "No" & use_data_loaded == "Yes"){
      vars_option <<- "READ"
      
      annot_vars <- colnames(all_read_part$mat_variables.0)
      annot_vars <- annot_vars[-which(annot_vars=="ID")]
      
      sel_list_1 <- list("Annotation variables" = annot_vars)
      sel_list_2 <- list("Annotation variables" = annot_vars)
      
      output$create_model_param <- renderUI({
        div(
          HTML("<h4>Choose your variables:</h4>"),
          
          fluidRow(
            column(6,
                   #variable to define sample groups in model
                   selectizeInput("var_def_groups", "Variable to DEFINE groups", choices= sel_list_1, multiple=FALSE, options = list(maxOptions=3000))
                   
            ),
            column(6,
                   #variable (or variables...) to classify sample groups defined by 'var_def_groups'
                   selectizeInput("var_to_classify", "Variable/s to CLASSIFY groups", choices= sel_list_2, multiple=TRUE, selected= sel_list_2[[1]][1], options = list(maxOptions=3000))
            )
          )
        )
      })
      
      subs_variable_group <- eventReactive(input$var_def_groups, {
        group_variable <- as.character(input$var_def_groups)
        class_var <- as.character(all_read_part$variables_info_mat[which(all_read_part$variables_info_mat[,1]==group_variable),2])
        if (class_var=="categoric") {
          legend_var <- sort(unique(as.character(all_read_part$mat_variables.0[,colnames(all_read_part$mat_variables.0)==group_variable])))
          div(
            HTML(paste("<i>",length(legend_var), " groups defined:</i>", sep="")),
            checkboxGroupInput("sub_group_var", "", choices = legend_var, selected = legend_var)
          )
        } else if (class_var=="numeric"){
          values <- round(as.numeric(as.character(all_read_part$mat_variables.0[,colnames(all_read_part$mat_variables.0)==group_variable])),3)
          div(
            HTML(paste("<i>",length(unique(values)), " groups defined:</i>", sep="")),
            sliderInput("slider_group_var_model", "", min=min(values, na.rm = T), max=max(values, na.rm = T), value = c(min(values, na.rm = T),max(values, na.rm = T))),
            sliderInput("slider_group_var_model_2", "", min=min(values, na.rm = T), max=max(values, na.rm = T), value = c(min(values, na.rm = T),max(values, na.rm = T))) 
          )
        }
        
      })
      
      
      output$create_model_param_2 <- renderUI({
        fluidRow(
          column(6,
                 subs_variable_group(),
                 div(style="text-align:left",actionBttn("button_create_model", "Create new model", style="fill", size="md", color="warning")),
                 busyIndicator(text="Running", wait=200)
          ),
          column(6
          )
        )
        
      })
      
    } 
    # else {
    #   output$create_model_param <- renderUI({
    #     div(
    #       HTML("<i style=color:red> Come on.. You have to select something!!</i>")
    #     )
    #   })
    # }
    
    
  })
  ########### MY MODEL param##########
  output$my_model_param <- renderUI({
    div(
      fluidRow(
        column(6,
               fileInput("my_model_browse", "Load your model", width="100%",
                         accept = c(
                           ".RData"),
                         buttonLabel = "Browse",
                         placeholder = ".RData")
        )
      )
      
      
    )
    
    
  })
  ############################
  ########WHEN BUTTON CREATE MODEL#########
  observeEvent(input$button_create_model, {
    if (input$var_def_groups == input$var_to_classify){
      output$error_model <- renderUI({
        HTML("<i style=color:red> Variable 'to CLASSIFY' can not be the same one as 'to DEFINE'!</i>")
      })
    } else {
      
      
      #####
      ##### CLASSIFICATORY VARS
      classy_vars <- input$var_to_classify
      classificatory_vars <- c()
      if (length( grep("All regions",classy_vars) )==1) {
        
        if (length( grep("Specific region/s",classy_vars) )==1){
          output$error_model <- renderUI({
            HTML("<i style=color:gray> Selecting 'All regions' and not considering 'Specific region/s'...</i>")
          })
          classy_vars <- classy_vars[-grep("Specific region/s",classy_vars)]
        }
        
        classificatory_vars <- as.character(all_profile_part[["segmented_file"]][,4])
        
        classy_vars <- classy_vars[-grep("All regions",classy_vars)]
        
      }
      
      if (length( grep("Specific region/s",classy_vars) )==1) {
        classificatory_vars <- c(classificatory_vars, as.character(input$spec_regions))
        
        classy_vars <- classy_vars[-grep("Specific region/s",classy_vars)]
        
      } 
      
      classificatory_vars <- c(classificatory_vars, classy_vars)
      
      #output$classy_vars <- renderUI(HTML(paste(classificatory_vars,sep="")))
      
      #####
      ##### GROUP BY VAR
      groups_by_var <- as.character(input$var_def_groups)
      values_var <- GLOBAL_DF[,which(colnames(GLOBAL_DF)==groups_by_var)]
      class_var <- class(values_var)
      
      # if (class_var=="numeric" | class_var=="integer") {
      if (!is.null(input$slider_group_var_model)) {
        values1 <- as.numeric(as.character(input$slider_group_var_model))
        min1 <- min(values1, na.rm = T)
        max1 <- max(values1, na.rm = T)
        values_var <- as.numeric(as.character(values_var))
        rows_in_annotation1 <- which(values_var>=min1 & values_var<=max1)
        
        values2 <- as.numeric(as.character(input$slider_group_var_model_2))
        min2 <- min(values2, na.rm = T)
        max2 <- max(values2, na.rm = T)
        rows_in_annotation2 <- which(values_var>=min2 & values_var<=max2)
        
        rows_in_annotation <- unique(c(rows_in_annotation1, rows_in_annotation2))
        values_var <- sort(values_var[rows_in_annotation])
        groups_var <- unique(values_var)
      }
      # if (class_var=="categoric" | class_var=="factor"){
      if (!is.null(input$sub_group_var)){
        groups_var <- input$sub_group_var
        rows_in_annotation <- which(GLOBAL_DF[,groups_by_var]%in%groups_var)
        values_var <- as.character(GLOBAL_DF[rows_in_annotation,groups_by_var])
      }
      n.groups <- length(groups_var)
      
      if (n.groups < 2) {
        output$error_model <- renderUI({
          HTML("<i style=color:red> Need at least two groups to do classification!</i>")
        })
      } else {
        
        output$error_model <- renderUI({
          HTML("")
        })
        
        if (class_var == "integer" | class_var == "numeric"){
          values_var <- as.numeric(as.character(values_var))
        }
        
        samples_in_groups <- GLOBAL_DF[rows_in_annotation,which(colnames(GLOBAL_DF)=="ID")[1]]
        
        mat <- as.data.frame(GLOBAL_DF[rows_in_annotation,classificatory_vars])
        rownames(mat) <- GLOBAL_DF[rows_in_annotation,which(colnames(GLOBAL_DF)=="ID")[1]]
        colnames(mat) <- classificatory_vars
        
        
        #########################################
        ######### RANDOM-FOREST MODEL ###########
        #########################################
        
        
        
        data_2 <- t(mat)
        
        
        fun_isNa_in_samples <- function(x){
          l_na <- length( which(is.na(x)) )
          if (l_na > 0) {
            ans <- TRUE
          } else {
            ans <- FALSE
          }
          ans
        }
        
        NAs_in_samples <- which(apply(data_2, 2, fun_isNa_in_samples)==TRUE)
        
        if ( length( NAs_in_samples ) > 0 ) {
          samples_in_groups <- samples_in_groups[-NAs_in_samples]
          data_2 <- data_2[,which(colnames(data_2)%in%samples_in_groups)]
        }
        
        
        rownames(data_2) <-paste("V",1:nrow(data_2),sep="")
        colnames(data_2) <- paste("S",1:ncol(data_2),sep="")
        data_2 <- rbind(data_2, rep(NA, length=dim(data_2)[2])) # in order to work when 'classificatory_vars' is just one...
        nrow_data <- nrow(data_2)
        groups_var <- sort(groups_var)
        
        if (length(which(is.na(values_var)))>0){
          samples_in_groups <- samples_in_groups[-which(is.na(values_var))]
          data_2 <- data_2[,-which(is.na(values_var))]
          values_var <- values_var[-which(is.na(values_var))]
          groups_var <- sort(unique(as.character(values_var)))
          n.groups <- length(groups_var)
        }
        
        
        
        
        
        s.part.pro <-  round(length(samples_in_groups)/n.groups)# total de samples dividit per el nombre de grups = s.part.pro
        min.count.groups <- min(table(as.factor(as.character(values_var)))) # n samples del grup amb menys samples
        
        
        
        
        
        if (min.count.groups < s.part.pro){
          s.part.pro <- round(min.count.groups * 0.75) # 75% de N samples del grup més petit
        } else if (min.count.groups == (length(samples_in_groups))/n.groups){
          s.part.pro <- round(min.count.groups * 0.75) # 75% de N samples del grup més petit
        }
        
        #### All levels included in data?
        
        e_vector <- c()
        e_class <- c()
        for (jj in 1:ncol(mat) ){
          name <- colnames(mat)
          all_levels_classy <- unique(mat[,jj])
          class_class <- class(mat[,jj])
          e_class <- c(e_class, class_class)
          
          if (length( all_levels_classy) > length(groups_var) & class_class != "numeric") {
            e_vector <- c(e_vector, "too_many")
            
          } else if (class_class != "numeric"){
            e_vector <- c(e_vector, "ok")
          }
        }
        
        
        
        
        # if ( (s.part.pro < (length(samples_in_groups) * 0.001)) |  ((length(samples_in_groups) * 0.001) < 1) ) {
        if ( (s.part.pro <= 1) |  (length(samples_in_groups) < 20) ) {
          
          output$error_model <- renderUI({
            HTML("<i style=color:red>Can't perform model!! N samples between groups have to be comparable. (Try other variables)</i>")
          })
          
          message <- print("Can't perform model!!")
          
          output$model_results_MAIN <- renderUI({HTML("")})
          output$vars_in_model <- renderUI({HTML("")})
          output$model_accuracy <- renderUI({HTML("")})
          output$final_stats_model_MAIN <- renderUI({HTML("")})
          output$final_stats_model <- renderTable({})
          output$n_samples_table <- renderTable({})
          output$bar_plot_predictions_MAIN <- renderUI({HTML("")})
          output$bar_plot_predictions <- renderPlotly({ plotly_empty() })
          output$bar_plot_real_MAIN <- renderUI({HTML("")})
          output$bar_plot_real <- renderPlotly({ plotly_empty() })
          output$final_table_model <- renderDataTable({})
          
        } else if (length(which(e_vector == "too_many"))>0) {
          
          output$error_model <- renderUI({
            HTML("Too many levels!<br> Try create a model without it... <br>")
          })
          
          message <- print("Can't perform model!!")
          
          output$model_results_MAIN <- renderUI({HTML("")})
          output$vars_in_model <- renderUI({HTML("")})
          output$model_accuracy <- renderUI({HTML("")})
          output$final_stats_model_MAIN <- renderUI({HTML("")})
          output$final_stats_model <- renderTable({})
          output$n_samples_table <- renderTable({})
          output$bar_plot_predictions_MAIN <- renderUI({HTML("")})
          output$bar_plot_predictions <- renderPlotly({ plotly_empty() })
          output$bar_plot_real_MAIN <- renderUI({HTML("")})
          output$bar_plot_real <- renderPlotly({ plotly_empty() })
          output$final_table_model <- renderDataTable({})
          
        } else {
          output$error_model <- renderUI({
            HTML("")
          })
          
          showElement("div_results_model")
          
          
          acc<-vector()
          
          iters <- 50
          ta <- 0
          
          list_model_preds <- list()
          list_ta <- list()
          for (i in 1:iters) {
            
            e_list <- list()
            for (ii in 1:n.groups){
              n.pro <- groups_var[ii]
              if (is.na(n.pro)){
                n.pro <- "na"
              }
              items.pro <- which(values_var==n.pro)
              l.pro <- length(items.pro)
              sub.pro <- sample(items.pro,s.part.pro)
              
              e_list[[paste(n.pro)]] <- sub.pro
            }
            
            training <- as.vector(do.call(cbind, e_list))
            
            data <- as.matrix(t(data_2[,training])[,-nrow_data])
            model<-randomForest(factor(values_var)[training] ~ . , data)
            
            pred<-predict(model, newdata=as.matrix(t(data_2[,-training])[,-nrow_data]))
            
            tabl<-table(pred, values_var[-training])
            tabl <- tabl[sort(rownames(tabl)), sort(colnames(tabl))]
            
            list_ta[[i]] <- tabl
            ta <- ta + tabl
            acc[i]<-sum(diag(tabl))/sum(tabl)
            list_model_preds[[i]] <- model$predicted
            print(i)
          }
          
          
          
          
          ta
          sum(diag(ta))/sum(ta)
          
          sen.spe <- matrix(0,ncol=n.groups,nrow=2)
          rownames(sen.spe)<-c("Sensitivity", "Specificity")
          colnames(ta)<-groups_var
          
          for (i in 1:ncol(ta)) {
            name <- colnames(ta)[i]
            m <- matrix(NA, ncol=2, nrow=2)
            colnames(m) <- c(name, "No Event")
            rownames(m) <- c(name, "No Event")
            m[1,1] <- ta[i,i]
            m[2,1] <- sum(ta[-i,i])
            m[1,2] <- sum(ta[i,-i])
            m[2,2] <- sum(ta)-(m[1,1] + m[2,1] + m[1,2])
            print(m)
            
            sen.spe[1,i] <- sensitivity(as.table(m))
            sen.spe[2,i] <- specificity(as.table(m))
            print(sen.spe)
            
            #                Reference
            # Predicted    Event    No Event
            #    Event      A         B
            # No Event      C         D
            #
            # The formulas used here are:
            #
            #              Sensitivity = A/(A+C)
            #
            #              Specificity =D/(B+D)
          }
          
          
          ta <- round(rbind(ta,sen.spe),3)
          mean_acc <- c(mean(acc),rep(NA, (dim(ta)[2]-1)))
          
          ans <- rbind(ta, mean_acc)
          ans
          
          #	NAME <- paste(out_dir, "Random_Forest_table_", name_SEG_FILE, ".txt", sep="")
          #	write.table(ans, NAME, sep="\t", col.names=NA, row.names=T, quote=F)
          
          ## Calcular percentatge de predicció del grup
          pre_mat <- matrix(NA,nrow=dim(data_2)[2],ncol=iters)
          rownames(pre_mat)<-colnames(data_2)
          colnames(pre_mat)<-1:iters
          
          for (ii in 1:length(list_model_preds)){
            x <- list_model_preds[[ii]]
            
            for (pp in 1:length(x)){
              n_x <- as.character(x[[pp]])
              s_x <- names(x)[pp]
              
              pre_mat[which(rownames(pre_mat)==s_x),ii] <- paste(n_x)
            }
            
          }
          
          list_group_counts <- apply(pre_mat,1,table)
          
          if (class(list_group_counts) == "matrix") {
            x <- list_group_counts
            list_group_counts <- split(x, rep(1:ncol(x), each = nrow(x)))
            for (jj in 1:length(list_group_counts)){
              names(list_group_counts[[jj]]) <- rownames(x)
            }
            names(list_group_counts) <- paste("V",1:length(list_group_counts),sep="")
          }
          
          confidence_th <- 0.75
          
          freq_mat <- as.data.frame(matrix(NA,nrow=dim(data_2)[2],ncol=n.groups+2))
          rownames(freq_mat)<-colnames(data_2)
          colnames(freq_mat)<- c(levels(factor(values_var)),"prediction_group", paste("confidence", "_", confidence_th, sep=""))
          
          freq_mat$sample_group <- values_var
          freq_mat <- cbind(samples_in_groups, freq_mat)
          colnames(freq_mat)[1] <- "sample_name"
          
          for (ii in 1:length(list_group_counts)){
            row <- list_group_counts[[ii]]
            sum_row <- sum(row)
            for (rr in 1:n.groups){
              group <- levels(factor(values_var))[rr]
              
              if (length(which(names(row)==group))==0){
                freq_mat[ii,which(colnames(freq_mat)==group)] <- 0
              } else if (length(which(names(row)==group))>0){
                freq <- as.numeric(row[which(names(row)==group)]/sum_row)
                
                freq_mat[ii,which(colnames(freq_mat)==group)] <- round(freq,3)
              }
            }
            
            freq_mat[ii,"prediction_group"] <- names(sort(freq_mat[ii,2:(2+n.groups-1)], decreasing=T))[1]
            
            if (freq_mat[ii,"prediction_group"]==freq_mat[ii, "sample_group"]){
              freq_mat$prediction[ii] <- "good"
            } else {
              freq_mat$prediction[ii] <- "bad"
            }
            
            
            first_pred <- sort(freq_mat[ii,2:(2+n.groups-1)], decreasing=T)[1]
            if (first_pred>=confidence_th){
              freq_mat[ii,grep("confidence", colnames(freq_mat))] <- "high"
            } else {
              freq_mat[ii,grep("confidence", colnames(freq_mat))] <- "low"
            }
          }
          
          rownames(data_2) <- c(classificatory_vars, "NA")
          freq_final <- cbind(freq_mat,t(data_2))
          freq_final <- freq_final[,-ncol(freq_final)]
          
          
          freq_tabl <- table(freq_final$prediction_group, freq_final$sample_group)
          freq_tabl
          sum(diag(freq_tabl))/sum(freq_tabl)
          
          
          ### functions
          v_by_group <- function(x, y, w, z){
            # x is item in list
            # y is matrix
            # w is column of list
            # z is column with values
            sub_y <- y[which(y[,w]==x),z]
            ans <- sub_y
            ans
          }
          fun.count <- function(x, y){
            # x is item in list
            # y is vector of names (paired with names(list))
            v_empty <- rep(NA, length(y))
            for (i in 1:length(y)){
              name_y <- y[i]
              n_name_y <- length(which(x==name_y))
              #print(n_name_y)
              v_empty[i] <- n_name_y
            }
            ans <- v_empty
            ans
          }
          
          #########################################
          ############# OUTPUTS Random Forest###############
          output$model_results_MAIN <- renderUI({
            HTML("<h3>RESULTS OF YOUR MODEL</h3>")
          })
          
          output$vars_in_model <- renderUI({
            HTML("Classifying <b>", groups_by_var, "</b> sample groups [<i>", paste(groups, collapse = ", "), "</i>] by <b>", as.character(input$var_to_classify), "</b> variable/s.")
          })
          # output$summary_data_model <- renderPrint(
          #   summary(data_2)
          # )
          
          output$model_accuracy <- renderUI({
            HTML("<table style='text-align:center'>
                 <tr>
                 <td><b>Global accuracy = </b></td>
                 <td>&nbsp&nbsp<b style='color:orange'>", round((ans[nrow(ans),1])*100,3), " (%)</b></td>
                 </tr>
                 </table>")
          })
          
          
          output$final_stats_model_MAIN <- renderUI({
            HTML("<h5 style=text-align:left> Classifier model stats:</h5>")
          })
          output$final_stats_model <- renderTable({
            ans1 <- ans * 100
            ans2 <- as.data.frame(cbind(rownames(ans1),ans1))
            colnames(ans2)[1] <- ""
            ans3 <- ans2[c("Sensitivity", "Specificity"),]
            ans3[,1] <- c("Sensitivity (%)", "Specificity (%)")
            
            ans3
          })
          
          
          
          data <- freq_final
          groups <- sort(groups_var)
          colors <- fixed_colors[1:length(groups)]
          
          mat <- mclapply(groups, v_by_group, y=data, w="sample_group",z="prediction_group")
          names(mat) <- groups
          mat2 <- as.data.frame(do.call(cbind, mclapply(mat, fun.count, y=groups)))
          rownames(mat2) <- groups
          # mat2 == cols are REAL groups / rows are PREDICTED groups
          
          total_samples_by_groups_real <- apply(mat2,2,sum)
          total_samples_by_groups_predicted <- apply(mat2,1,sum)
          
          mat_samples <- rbind(total_samples_by_groups_real, total_samples_by_groups_predicted)
          rownames(mat_samples) <- c("Real samples (n)", "Predicted samples (n)")
          
          total_samples <- sum(total_samples_by_groups_real)
          
          
          # output$n_samples_MAIN <- renderUI({
          #   HTML("<h5 style=text-align:left> N samples by groups:</h5>")
          # })
          output$n_samples_table <- renderTable(rownames=T,{
            mat_samples
          })
          
          #### Accuracy bar
          # mt <- as.data.frame(t(c("", round((rev(as.vector(table(data[,"prediction"])))/total_samples)*100,2)  )))
          # colnames(mt) <- c("acc", "good", "bad")
          #
          # bar_acc <- plot_ly(y=as.character(mt[,"acc"]), x=as.numeric(as.character(mt[,"good"])), name="good",type = "bar", marker=list(color="green"), showlegend=F, orientation='h')%>%
          #   add_trace(x=as.numeric(as.character(mt[,"bad"])), name="bad", marker=list(color="red"), orientation='h')%>%
          #   layout(
          #     barmode="stack",
          #     xaxis=list(visible=F),
          #     yaxis=list(visible=F, title=c("%"), titlefont=list(size=30)),
          #     font = list(size=2),
          #     margin = list(l = 10, r = 5, t = 0, b = 0),
          #     height=75
          #   )
          # output$acc_bar <- renderPlotly({
          #   bar_acc
          # })
          
          
          ###### Bar plot predictions####
          mat3 <- round((mat2/total_samples_by_groups_predicted)*100,0)
          counts <- as.data.frame(cbind(groups, mat3))
          colnames(counts)[1] <- c("groups")
          
          bar_pred <- plot_ly(counts, x=counts[,"groups"], y=counts[,groups[1]], name=groups[1],type = "bar", marker=list(color=colors[1]), showlegend=T)%>%
            layout(
              yaxis =list(title= "", type="linear", range=NULL),
              barmode="stack",
              
              font = list(size=20),
              margin = list(l = 90, r = 20, t = 20, b = 40)
            )
          
          for(i in 2:length(groups)){
            bar_pred <- add_trace(bar_pred, y=counts[,groups[i]], name=groups[i], marker=list(color=colors[i]), showlegend=T)
          }
          
          output$bar_plot_predictions_MAIN <- renderUI({
            HTML("<h4 style=text-align:center>Predicted-tagged samples between <u>real</u> <b><i>", groups_by_var, "</i></b> groups</h4>")
          })
          output$bar_plot_predictions <- renderPlotly({
            bar_pred
          })
          
          ###### Bar plot: real samples in predicted groups####
          mat3 <- round((t(mat2)/total_samples_by_groups_real)*100,0) # mat2 is transposed!!! and it is divided by total_samples_by_groups_real
          counts <- as.data.frame(cbind(groups, mat3))
          colnames(counts)[1] <- c("groups")
          
          bar_real <- plot_ly(counts, x=counts[,"groups"], y=counts[,groups[1]], name=groups[1],type = "bar", marker=list(color=colors[1]), showlegend=T)%>%
            layout(
              yaxis =list(title= "", type="linear", range=NULL),
              barmode="stack",
              
              font = list(size=20),
              margin = list(l = 90, r = 20, t = 20, b = 40)
            )
          
          for(i in 2:length(groups)){
            bar_real <- add_trace(bar_real, y=counts[,groups[i]], name=groups[i], marker=list(color=colors[i]), showlegend=T)
          }
          
          output$bar_plot_real_MAIN <- renderUI({
            HTML("<h4 style=text-align:center>Real-tagged samples between <u>predicted</u> <b><i>", groups_by_var, "</i></b> groups</h4>")
          })
          output$bar_plot_real <- renderPlotly({
            bar_real
          })
          
          output$final_table_model <- renderDataTable({
            freq_final
          },
          options = list(scrollX=TRUE, pageLength=15, lengthMenu=c(5,10,15,20), paging=TRUE, searching=TRUE, info=FALSE)
          )
          
          output$dw_button_freq_final <- renderUI(
            div(style="text-align:right",
                downloadButton("TSV_class_model", label="TSV")
            )
          )
          
          output$TSV_class_model <- downloadHandler(
            filename = paste("Classifier_model_predictions_", groups_by_var, "_by_", as.character(input$var_to_classify), "_", Sys.Date(), ".tsv", sep=""),
            content = function(file) {
              write.table(freq_final, file, sep="\t", col.names=T, row.names=F, quote=F, dec=".")
            }
          )
          
          #########################################
        }
        
      } #else (CAN PERFORM MODEL!)
      
      
      
      
      
      
    } # else of (input$var_to_classify == input$var_def_groups)
    
  }) #observeEvent(input$button_create_model)
  
  # observeEvent( input$button_param_create_model, {
  #   
  # })
  
  
  # ########WHEN BUTTON MY MODEL#########
  # observeEvent(input$button_my_model, {
  #   
  #   # hide("side_menu_model")
  #   # #hide("side_menu_create_model")
  #   # output$side_menu_my_model <- renderMenu(
  #   #   sidebarMenu(
  #   #     menuItem(HTML("<b style=color:#00cc66>USE MY MODEL</b>"), startExpanded =T,
  #   #              menuSubItem(HTML("Change mode to <b style=color:orange> &nbsp CREATE </b>"), tabName="home_model", icon=icon("microchip")),
  #   #              menuSubItem("Parameters", tabName="tab_my_model_param", icon=icon("cogs"), selected = T),
  #   #              menuSubItem("Results modeling", tabName="tab_my_model_results", icon=icon("object-ungroup"))
  #   #     )
  #   #   )
  #   # )
  #   # #show("side_menu_my_model")
  #   
  #   
  #   
  #   
  #   
  # })
  
  
  
  
  
  
  
  
  
  
} #function (input,output)


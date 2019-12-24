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
library(Cairo)
library(DT)
# setwd("/srv/shiny-server/CNApp/")


dashboardPage(title="CNApp - Copy Number Alterations Integrative Analysis",
  #theme = shinytheme("paper"),
  #shinythemes::themeSelector(), # easy select shiny themes (comment when selection is set!!!)
  dashboardHeader(

    title=HTML("<h2 style='text-align:left;margin:5px 0 0 20px'>CNApp</h2> <img style='margin:-60px 0 0 50px' src='ITEMS/imatges/chrom_cnapp.png' width='40%' heigth='40%'>"),
                  
                  dropdownMenu(type="messages", icon=icon("home"), headerText = "", badgeStatus = NULL,
                               messageItem(icon=icon("home"), href="https://tools.idibaps.org/CNApp/home.html",
                                           from= tags$div(tags$h4("CNApp web page")),
                                           message=HTML("")
                               )
                  ),
                  
                  dropdownMenu(type="messages", icon=icon("info-circle"), headerText = "", badgeStatus = NULL,
                               messageItem(icon=icon("info-circle"), href="https://tools.idibaps.org/CNApp/#about",
                                           from= tags$div(tags$h4("About CNApp")),
                                           message=HTML("")
                               )
                  ),
                  
                  dropdownMenu(type="messages", icon=icon("github"), headerText = "", badgeStatus = NULL,
                               messageItem(icon=icon("github"), href="https://github.com/ait5/CNApp",
                                           from= tags$div(tags$h4("github.com/ait5/CNApp")),
                                           message=HTML("")
                               )
                  ),
                  
                  dropdownMenu(type="messages", icon=icon("refresh"), headerText = "CNApp will refresh and lose your session...", badgeStatus = NULL,
                               messageItem(icon=icon("exclamation-triangle"), href=NULL,
                                           from= tags$div(tags$h4("Are you sure?"),br()),
                                           message=actionBttn("clear", "Yes, clear all!", style="simple", size="sm", color="primary")
                               )
                  )

  ),
  
  ######dashboardSidebar#######
  dashboardSidebar(div(id="dashboard_sidebar", style="position:fixed; margin: 0 5px 0 5px; width: 220px; min-width: 200px; background-color:#253744;",
    div(id="body_sidebar",
    sidebarMenu(
      menuItem(text = div(HTML("<i class='fa fa-desktop' style='color:white'></i>&nbsp&nbsp HOME &nbsp&nbsp&nbsp&nbsp&nbsp&nbsp")), tabName="wellcome_tab"),
      menuItem(text = div(HTML("<i class='fa fa-folder-open' style='color:white'></i>&nbsp&nbsp LOAD DATA &nbsp&nbsp&nbsp&nbsp&nbsp&nbsp")), tabName="load_data")
    ),
    sidebarMenuOutput("side_menu_cin"),
    sidebarMenuOutput("side_menu_cin.1"),
    sidebarMenuOutput("side_menu_profiling"),
    sidebarMenuOutput("side_menu_profiling.1"),
    sidebarMenuOutput("side_menu_model")
    # div(style = "position: fixed; overflow: visible; padding-top:15em;", sidebarMenuOutput("side_menu_create_model")),
    # div(style = "position: fixed; overflow: visible; padding-top:15em;", sidebarMenuOutput("side_menu_my_model"))
    

    ),

    div(id="footer", style="margin-left: 0px; margin-bottom: 0px; text-align: center; position: absolute; top:425px; background-color:#253744; padding: 50px 0 20px 0;",
    HTML("
         
         
							<a href='http://www.idibaps.org/' class='icons'><img src='ITEMS/imatges/sidebar/IDIBAPS1.png', alt='idibaps-image-link', width='150px', height='50px'></img></a>
<br><br><br>
							<a href='http://www.ciberehd.org/' class='icons'><img src='ITEMS/imatges/sidebar/ciberehd.png', alt='ciberehd-image-link', width='150px', height='50px'></img></a>
<br><br><br>
						<p class='copyright'>&copy; CNApp web application. Designed by <a href=''>Sebastià Franch</a></p>.
					
         
         ")
    )
    #includeHTML("www/en/sidebar_page.html")

  )),
  
  ##########dashboardBody#############
  dashboardBody(
    HTML("<script async src='https://www.googletagmanager.com/gtag/js?id=UA-121173624-1'></script><script>window.dataLayer = window.dataLayer || [];  function gtag(){dataLayer.push(arguments);}  gtag('js', new Date());  gtag('config', 'UA-121173624-1');</script>"
    ),
    useShinyjs(),
    
    tabItems(
      ######Tab WELCOME PAGE######
      tabItem(tabName="wellcome_tab",
              
             includeHTML("www/en/wellcome_page.html")
              
              
      ), #tabName="load_data"
      
      ######Tab LOAD DATA######
      tabItem(tabName="load_data",

              fluidRow(
                column(4,fileInput("data_browse", "Load your data", width="100%",
                                   accept = c(
                                     "text/csv",
                                     "text/comma-separated-values,text/plain",
                                     ".csv",
                                     ".bed",
                                     ".tsv",
                                     ".txt/tab-separeted-values"),
                                   buttonLabel = "Browse")
                ),
                column(1,
                       div(style="position:relative; top:2.5em; text-align:left;", title="empty data", actionLink("empty", label=icon("trash", "fa-1x")))
                ),
                column(2,
                       br(),
                       div(style="text-align:left;", checkboxInput("use_demo", "Demo data", value=F))
                       ),
                column(1),
                column(2, 
                       #Help menu for format of input file
                       div(style="text-align:right;", actionLink("helpformat","Help", icon=icon("question-circle-o"))),
                       bsModal("modal","HELP: Input file","helpformat", includeHTML("./aux_files/help_with_input.html"))
                       

                )
              ),
              
                  fluidRow(
                    column(9,
                           
                           # File specifications panel
                           uiOutput("file_specifics"),
                           
                           textOutput('error_demo_read'),
                           tags$head(tags$style("#error_demo_read{color: red;
                                                  font-size: 10px;
                                                  font-style: italic;
                                                  }"
                           )
                           ),
                           
                           textOutput('error_reading_data'),
			   uiOutput('error_too_many_samples'),
                           textOutput('error_reading_annot_data'),
                           
                           tags$head(tags$style("#error_reading_data{color: red;
                                      font-size: 10px;
                                      font-style: italic;
                                      }"
                           )
                           ),
                           br(),
                           uiOutput("go_to_funs")
                    )
                  ),
                  br(),
                  hr(),

                  uiOutput("summary_data_MAIN"),
                  uiOutput("dim_data"),
                  verbatimTextOutput("summary_data")



              
      ), #tabName="load_data"
      
      ########Tab PARAM CIN########
      tabItem(tabName="user_param_cin",
              HTML("<h3><i class='fa fa-star' style='color:orange'></i>&nbsp&nbsp RE-SEG & SCORE</h3>"),
              br(),
              
              div(id="USER_PARAM_cin",
                  div(style="text-align:right;", uiOutput("help_param_cin")),
                  uiOutput("user_parameters_cin"),
                  uiOutput("user_parameters_cin_2"),
                  uiOutput("user_parameters_skip_reseg"),
                  uiOutput("user_parameters_cin_3")
              ),
              
              uiOutput("anal_run_cin")
              
      ),#tabName="user_param_cin"
      
      ########Tab RESULTS CIN#######
      tabItem(tabName="results_cin",
              HTML("<h3><i class='fa fa-star' style='color:orange'></i>&nbsp&nbsp RE-SEG & SCORE</h3>"),
              br(),
                        
              tabsetPanel(id="tab",
                          
                          #Plots
                          tabPanel("Re-segmentation",

                                   br(),
                                   br(),

                                   HTML("<h4 style='text-align:center'><b>Download re-segmented samples file:</b></h4>"),
                                   div(style="text-align:center", downloadButton("ann_segments_all", label="TSV")),
                                   div(style="text-align:center", uiOutput("samples_down_in_re_seg")),
                                   
                                   br(),
                                   br(),
                                   br(),
                                   
                                   #INDIVIDUAL PLOT:
                                   HTML("<h4 style='text-align:center'><b>Individual segments plots</b></h4>"),
                                   HTML("<h5 style='text-align:center'>(<i style=color:black>before</i> / <i style=color:green>after</i> re-segmentation)</h5>"),
                                   uiOutput("choose_sample_to_plot"),
                                   
                                   plotOutput("bef_before", height = "250px"),
                                   plotOutput("bef_after", height = "250px"),
                                   uiOutput("dw_buttons_plots_indv_segs"),
                                   
                                   br(),
                                   br(),
                                   br(),
                                   
                                   #Freq CN plots:
                                   HTML("<h4 style='text-align:center'><b>Recurrent CNAs & frequency plot</b></h4>"),
                                   div(style="text-align:center", uiOutput("cnfreq_options")),
                                   fluidRow(
                                     column(2),
                                     column(4,
                                            uiOutput("sub_variable")
                                     ),
                                     column(4,
                                            uiOutput("sub_sub_variable")
                                     ),
                                     column(2)
                                   ),
                                   
                                   uiOutput("cnfreq_plot_button"),
                                   br(),
                                   uiOutput("few_samples"),
                                   plotlyOutput("freq_cn_filt"),
                                   uiOutput("dw_button_cnfreq"),
                                   br(),
                                   br()

                                   
                          ),
                          
                          
                          tabPanel("CNA Scores",value="scores",
                                   br(),
                                   br(),
                                  
                                   #Boxplot CNA Scores:
                                   HTML("<h4 style='text-align:center'><b>CNA Scores distribution</b></h4>"),
                                   div(style="text-align:center", uiOutput("box_scores_options")),
                                   fluidRow(
                                     column(2),
                                     column(4,
                                            uiOutput("sub_variable_score")
                                     ),
                                     column(4,
                                            uiOutput("sub_sub_variable_score")
                                     ),
                                     column(2)
                                   ),
                                   
                                   uiOutput("score_plot_button"),
                                   br(),
                                   # uiOutput("few_samples"),
                                   
                                   # fluidRow(
                                   #   column(4,
                                   #          plotOutput("box_BCS", height="300px")
                                   #          # uiOutput("dw_button_BCS")
                                   #   ),
                                   #   column(4,
                                   #          plotOutput("box_FCS", height="300px")
                                   #   ),
                                   #   column(4,
                                   #          plotOutput("box_GCS", height="300px")
                                   #   )
                                   # ),
                                   fluidRow(
                                     column(2),
                                     column(8,
                                            plotOutput("box_BCS"),
                                            uiOutput("dw_button_BCS_plot")
                                            ),
                                     column(2)
                                   ),
                                   br(),
                                   
                                   fluidRow(
                                     column(2),
                                     column(8,
                                            plotOutput("box_FCS"),
                                            uiOutput("dw_button_FCS_plot")
                                     ),
                                     column(2)
                                   ),
                                   br(),
                                   
                                   fluidRow(
                                     column(2),
                                     column(8,
                                            plotOutput("box_GCS"),
                                            uiOutput("dw_button_GCS_plot")
                                     ),
                                     column(2)
                                   ),
                                  
                                   
                                   br(),
                                   br(),
                                   br(),
                                   
                                   HTML("<h4 style='text-align:center'><b>CNA Scores and annotations by samples</b></h4>"),
                                   br(),
                                   fluidRow(
                                     column(11,
                                            dataTableOutput("scores"),
                                            div(style="text-align:right", downloadButton("download_scores",label="TSV"))
                                            )
                                   ),

                                   br(),
                                   # br(),
                                   # br(),
                                   # 
                                   # 
                                   # 
                                   # 
                                   # HTML("<h4 style='text-align:center'><b>Hierarchical clustering for sample CNA Scores</b></h4>"),
                                   # uiOutput("ht_scores_plot_button"),
                                   # br(),
                                   # plotlyOutput("heatmap_scores",width="100%", height="600px"),
                                   # uiOutput("dw_buttons_ht_scores"),
                                   # 
                                   # br(),
                                   br()
                                   # downloadButton("download_scores_raw",label="Download scores (raw)")
                          ),
                          
                          tabPanel("Variable association",value="clinical",
                                   
                                   br(),
                                   HTML("<h4><b>Parametric tests</b> <i>(p-values)</i></h4>"),
                                   dataTableOutput("pval.parametric"),
                                   div(style="text-align:right", downloadButton("download_pval.parametric",label="TSV")),
                                   
                                   br(),
                                   br(),
                                   br(),
                                   
                                   HTML("<h4><b>Non-parametric tests</b> <i>(p-values)</i></h4>"),
                                   dataTableOutput("pval.nonparametric"),
                                   div(style="text-align:right", downloadButton("download_pval.nonparametric",label="TSV")), 
                                   
                                   br(),
                                   br()
                          ),
                          
                          tabPanel("Survival analysis",value="clinical",
                                   
                                   br(),
                                   HTML("<h4><b>Survival analysis</b> <i>(by annotation variables)</i></h4>"),
                                   br(),
                                   HTML("<h5>Choose your variables:</h5>"),
                                   br(),
                                   
                                   div(style="text-align:center", uiOutput("prepare_surv_analysis")),
                                   uiOutput("sub_g_var"),
                                   uiOutput("not_enough_groups_surv"),
                                   br(),
                                  
                                   uiOutput("button_run_surv"),
                                   
                                   br(),
                                   
                                   fluidRow(
                                     column(2),
                                     column(8,
                                            # uiOutput("message"),
                                            plotOutput("surv_curves_plot"),
                                            uiOutput("dw_button_survival_plot")
                                     ),
                                     column(2)
                                   ),
                                   
                                   br(),
                                   br(),
                                   br()
                          )

                          # tabPanel("Predictions",value="predict",
                          #          downloadButton("download_predictions",label="Download predictions"),
                          #          dataTableOutput("predictions.tab")
                          # )
              )
              
      ),#tabName="results_cin"
      
      ########Tab PARAM PROFILING#######
      tabItem(tabName="user_param_seg",
              HTML("<h3><i class='fa fa-map' style='color:orange'></i>&nbsp&nbsp REGION PROFILE</h3>"),
              br(),

                     # User analasys parameters panel
                     div(id="USER_PARAM",
                         div(style="text-align:right;", uiOutput("help_param")),
                         uiOutput("user_parameters"),
                         
                         br(),
                         uiOutput("user_parameters_2"),

                         uiOutput("user_parameters_3")
                     ),
                     uiOutput("anal_run")
              
              ),#tabName="user_param_seg"
      
      ########Tab RESULTS PROFILING#########
      tabItem(tabName="seg_cna_profile_tab",
              HTML("<h3><i class='fa fa-map' style='color:orange'></i>&nbsp&nbsp REGION PROFILE</h3>"),
              br(),
              uiOutput("cn_profiling_TITLE"),
              
              sidebarLayout(position="right",
                            sidebarPanel(width=3,
                                         uiOutput("user_friendly_plot_param_1"),
                                         uiOutput("user_friendly_plot_param_2"),
                                         uiOutput("user_friendly_plot_param_3"),
                                         uiOutput("message_text")
                            ),
                            mainPanel(width=9,
                                      
                                      tabsetPanel(
                                        
                                        tabPanel(title="CNA region profiles",
                                                 div(style="text-align:center", htmlOutput("cna_profile_heatmap_MAIN")),
                                                 br(),
                                                 plotlyOutput("cna_profile_heatmap", height="800px"),

                                                 fluidRow(
                                                   column(6,
                                                          uiOutput("dw_annotation_tracks_tab1")
                                                   ),
                                                   column(6,
                                                          uiOutput("dw_buttons_cna_profile_htmp")
                                                   )
                                                 ),
                                                 
                                                 # br(),
                                                 # br(),
                                                 # 
                                                 # uiOutput("cor_test_annot_reg"),
                                                 
                                                 br(),
                                                 br(),
                                                 br(),
                                                 
                                                 div(style="text-align:center", htmlOutput("cna_profile_dendro_MAIN")),
                                                 uiOutput("cna_prof_dendro_plot_button"),

                                                 # uiOutput("prova_1"),
                                                 # uiOutput("prova_2"),
                                                 # uiOutput("prova_3"),
                                                 # uiOutput("prova_4"),

                                                 br(),
                                                 
                                                 # uiOutput("dendro_type_box"),
                                                 plotlyOutput("cna_profile_dendrogram", height="700px", width = "100%"),
                                                 uiOutput("dw_buttons_dendro_profile"),
                                                 
                                                 br(),
                                                 br(),
                                                 br(),
                                                 
                                                 uiOutput("tab_genes_MAIN"),
                                                 uiOutput("select_reg_genes"),
                                                 
                                                 
                                                 dataTableOutput("tab_genes_reg"),
                                                 uiOutput("dw_tab_genes_bttn"), 
                                                 
                                                 br(),
                                                 br()
                                                 
                                        ),
                                        
                                        tabPanel(title="CNA region frequencies",
                                                 div(style="text-align:center", htmlOutput("cna_events_by_regs_MAIN")),
                                                 br(),
                                                 plotlyOutput("cna_events_by_regs", height="600px"),
                                                 uiOutput("dw_buttons_cna_events_by_regs"), 
                                                 
                                                 
                                                 br(),
                                                 br(),
                                                 br(),
                                                 
                                                 div(style="text-align:center", htmlOutput("cna_events_MAIN")),
                                                 br(),
                                                 plotlyOutput("cna_events", height="700px"),
                                                 
                                                 fluidRow(
                                                   column(6,
                                                          uiOutput("dw_annotation_tracks_tab3")
                                                   ),
                                                   column(6,
                                                          uiOutput("dw_buttons_cna_events_barplot")
                                                   )
                                                 ),
                                                 
                                                 br(),
                                                 br()
                                                 
                                        ),
                                        
                                        tabPanel(title="Correlation profiles",
                                                 div(style="text-align:center", htmlOutput("unsupervised_corr_heatmap_MAIN")),
                                                 br(),
                                                 plotlyOutput("unsupervised_corr_sample_heatmap", height = "700px"),

                                                 fluidRow(
                                                   column(6,
                                                          uiOutput("dw_annotation_tracks_tab2")
                                                   ),
                                                   column(6,
                                                          uiOutput("dw_buttons_unsup_corr_htmp")
                                                   )
                                                 ), 
                                                 
                                                 br(),
                                                 br(),
                                                 br(),
                                        
                                                 
                                                 div(style="text-align:center", htmlOutput("cor_samples_dendro_MAIN")),
                                                 uiOutput("dendro_type_box_cor_plot_button"),
                                                 
                                                 br(),
                                                 uiOutput("dendro_type_box_cor"),
                                                 plotlyOutput("cor_samples_dendrogram", height="700px", width="100%"),
                                                 uiOutput("dw_buttons_dendro_cor"),
                                                 
                                                 br(),
                                                 br()
                                                 
                                        )
                                        
                                      )
                                      
                            )
              )
              
      ),#tabName="corr_tab"
      
      #######Tab LIMMA#######
      tabItem(tabName="limma_tab",
              HTML("<h3><i class='fa fa-map' style='color:orange'></i>&nbsp&nbsp REGION PROFILE</h3>"),
              br(),
              uiOutput("descriptive_regions_TITLE"),
              
              sidebarLayout(position="right",
                sidebarPanel(width=3,

                    uiOutput("side_param_limma"),
                    uiOutput("side_param_limma_2"),
                    uiOutput("signif"),
                    uiOutput("text_limma")

                ),
                mainPanel(width=9,


                          tabsetPanel(
                            tabPanel(title="Student's t-test",
                              
                                     div(style="text-align:center", htmlOutput("limma_heatmap_MAIN")),
                                     withTags(i(uiOutput("not_differences_2"))),
                                     plotlyOutput("limma_heatmap", height="900px"),
                                     uiOutput("dw_buttons_limma_heatmap"),
                                     
                                     div(style="text-align:center", htmlOutput("reg_boxplot_MAIN")),
                                     #withTags(i(uiOutput("not_differences"))),
                                     plotOutput("reg_boxplot", height="500px"),
                                     uiOutput("dw_buttons_reg_boxplot"),
                                     
                                     br(),
                                     div(style="text-align:left", htmlOutput("t_genes_reg_MAIN")),
                                     dataTableOutput("t_genes_reg"),
                                     uiOutput("dw_tab_genes_in_reg"), 
                                     
                                     br(),
                                     br()       
                                     
                            ),
                            
                            tabPanel(title="Fisher's test",
                                     
                                     div(style="text-align:center", htmlOutput("fisher_heatmap_MAIN")),
                                     withTags(i(uiOutput("not_differences_2_fisher"))),
                                     plotlyOutput("fisher_heatmap", height="900px"),
                                     uiOutput("dw_buttons_fisher_heatmap"),
                                     
                                     div(style="text-align:center", htmlOutput("reg_stacked_MAIN")),
                                     #withTags(i(uiOutput("not_differences_fisher"))),
                                     plotlyOutput("reg_stacked", height="500px"),
                                     uiOutput("dw_buttons_reg_stacked"),
                                     
                                     br(),
                                     div(style="text-align:left", htmlOutput("t_genes_reg_MAIN_fisher")),
                                     dataTableOutput("t_genes_reg_fisher"),
                                     uiOutput("dw_tab_genes_in_reg_fisher"), 
                                     
                                     br(),
                                     br()  
                              
                            )
                          )
                          
                          
                          
                )
              )

              
      ),#tabName="limma_tab"
      
      
     
      # ########Tab HOME MODEL#########
      # tabItem(tabName="home_model",
      #         div(tags$h3("User classifier model ", tags$i("( by Random Forest)"))),
      #         
      #         div(id="home_model",
      #             div(style="text-align:right;", uiOutput("help_home_model")),
      #             uiOutput("home_model")
      #         )
      #         
      # ),#tabName="home_model"

      
      
      #########Tab CREATE MODEL PARAM#########
      tabItem(tabName="tab_create_param",
              HTML("<h3><i class='fa fa-microchip' style='color:orange'></i>&nbsp&nbsp CLASSIFIER MODEL</h3>"),
              br(),
              fluidRow(
                column(10,
                       HTML("<h4>Generate your classifier model <i>(by Random Forest)</i></h4>")
                       ),
                column(2,
                       div(style="text-align:right;", actionLink("helpparam_machine","Help", icon=icon("question-circle-o"))),
                       bsModal("modal_param_machine","HELP:  User parameters CLASSIFIER MODEL","helpparam_machine", includeHTML("./aux_files/help_with_param_model.html"))
                )
              ),

              
              br(),
              uiOutput("run_parts_in_create_model_cin"),
              uiOutput("run_parts_in_create_model_profile"),
              uiOutput("run_parts_in_create_model_data_loaded"),
              
              # br(),
              # tableOutput("global_df"),
              
              br(),
              br(),
              br(),
              
              uiOutput("create_model_param"),
              # uiOutput("subselect_groups_model"),
              uiOutput("create_model_param_2"),
              uiOutput("error_model"), 
              
              # br(),
              # uiOutput("classy_vars"),
              # tableOutput("global_df_2"),
              #uiOutput("jjjj"),
              
              br(),
              br(),
              br(),
              
              hidden(div(id="div_results_model",
              uiOutput("model_results_MAIN"),
              
              br(),
              wellPanel(
              uiOutput("vars_in_model")),
              # verbatimTextOutput("summary_data_model"),
              br(),
              
              fluidRow(
                column(4,
                       div(style="text-align:left", wellPanel(uiOutput("model_accuracy")))
                       # br(),
                       # div(style="text-align:left", wellPanel(plotlyOutput("acc_bar", inline = FALSE)))
                )
              ),

              br(),

              wellPanel(
                uiOutput("final_stats_model_MAIN"),
                div(style="text-align:left", tableOutput("final_stats_model")),
                br(),
                div(style="text-align:left", tableOutput("n_samples_table"))
              ),

              
              br(),
              br(),
              br(),
              
              uiOutput("bar_plot_predictions_MAIN"),
              div(style="text-align:center", wellPanel(plotlyOutput("bar_plot_predictions", width = "80%", inline = TRUE))),
              
              br(),
              br(),
              br(),
              
              uiOutput("bar_plot_real_MAIN"),
              div(style="text-align:center", wellPanel(plotlyOutput("bar_plot_real", width = "80%", inline = TRUE))),
              
              br(),
              br(),
              br(),
              
              dataTableOutput("final_table_model"),
              
              uiOutput("dw_button_freq_final"),
              
              br(),
              br()
              ))#div
              
              
      )#tabName="model_results"
      # ),#tabName="model_results"
      # 
      # #########Tab CREATE MODEL RESULTS#########
      # # tabItem(tabName="tab_create_results",
      # #         div(tags$h3("Generate your classifier model ", tags$i("( by Random Forest)")))
      # #         
      # # ),#tabName="model_results"
      # 
      # 
      # 
      # #########Tab MY MODEL PARAM#########
      # tabItem(tabName="tab_my_model_param",
      #         HTML("<h3><i class='fa fa-microchip' style='color:orange'></i>&nbsp&nbsp CLASSIFIER MODEL</h3>"),
      #         br(),
      #         HTML("<h4>User classifier model <i>(by Random Forest)</i></h4>"),
      #         
      #         br(),
      #         uiOutput("my_model_param"),
      #         br(),
      #         
      #         div(style="text-align:center",actionBttn("button_my_model", "Use my model", style="fill", size="md", color="success")),
      #                       busyIndicator(text="Running", wait=200), 
      #         
      #         br(),
      #         br()
      #         
      # ),#tabName="model_results"
      # 
      # #########Tab MY MODEL RESULTS#########
      # tabItem(tabName="tab_my_model_results",
      #         div(tags$h3("User classifier model results", tags$i("( by Random Forest)")))
      #         
      # )#tabName="model_results"
      
    )#tabItems
    #################
  )#dashboardBody
  
)#dashboardPage

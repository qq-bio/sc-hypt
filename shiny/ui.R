# ui.R

library(shiny)
library(DT)
library(shinyWidgets)

setwd("/xdisk/mliang1/jingh/Rshiny/Shiny_Hypertension")

# File listings
seurat_files <- list.files("www/seurat_object", full.names = FALSE)
deg_files    <- list.files("www/deg_list/deg_download",    full.names = FALSE)

# Named vectors: names → description text
seurat_descriptions <- c(
    "C57BL6_HYP.rds"  = "Processed data of hypothalamus from C57BL/6 mouse.",
    "C57BL6_LV.rds"   = "Processed data of Left ventricle from C57BL/6 mouse.",
    "C57BL6_LK.rds"   = "Processed data of Left Kidney from C57BL/6 mouse.",
    "SP_HYP.rds"      = "Processed data of Hypothalamus from WKY rat and SHR rat.",
    "SP_LK.rds"       = "Processed data of Left Kidney from WKY rat and SHR rat.",
    "SP_LV.rds"       = "Processed data of Left ventricle from WKY rat and SHR rat.",
    "SP_MCA.rds"      = "Processed data of middle cerebral artery from WKY rat and SHR rat.",
    "SP_MSA.rds"      = "Processed data of 3rd mesenteric artery from WKY rat and SHR rat.",
    "SS_HYP.rds"      = "Processed data of Hypothalamus from SD rat and SS rat.",
    "SS_LK.rds"       = "Processed data of Left Kidney from SD rat and SS rat.",
    "SS_LV.rds"       = "Processed data of Left ventricle from SD Rat and SS Rat.",
    "SS_MSA.rds"      = "Processed data of 3rd mesenteric artery from SD Rat and SS Rat."
)

deg_descriptions <- c(
    "DEG.all.out"  = "DEG between control and treatment by wilcoxon-test.",
    "pseudo.DEG.all.out"  = "DEG between control and treatment by DESeq2 on pseudobulk data."
)

# UI definition
ui <- tagList(
    tags$head(
        # FontAwesome for icons
        tags$link(
            rel  = "stylesheet",
            href = "https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.5.0/css/all.min.css"
        ),
        # Global styles
        tags$style(HTML("
            
            /* Navbar styling */
            .navbar {
                position: fixed; top: 0; left: 0; right: 0;
                min-height: 70px; z-index: 9999;
                box-shadow: 0 2px 4px rgba(0,0,0,0.1);
            }
            .navbar-default {
                background-color: #b3c3d5;
                border-color: #b3c3d5;
            }
            .navbar-default .navbar-brand,
            .navbar-default .navbar-nav > li > a {
                color: #193D5B !important;
                font-size: 18px !important;
                font-weight: bold;
                display: flex;
                align-items: center;
                justify-content: center;
                height: 70px;
                gap: 6px;
            }
            .navbar-default .navbar-brand:hover,
            .navbar-default .navbar-nav > li > a:hover {
                background-color: #0056b3 !important;
                color: #ffffff !important;
            }
            .navbar .navbar-brand { float: left; margin-left: 20px; }
            .navbar .navbar-nav  { float: right; margin-right: 20px; }
            
            /* Body and container styling */
            body {
                overflow-x: hidden;
                background-color: #f5f5f5;
                padding-top: 70px;
            }

            .container-main,
            .container-main2 {
                max-width: 1400px;
                margin: 40px auto;
                background-color: white;
                padding: 20px 50px;
                box-shadow: 0 2px 4px rgba(0,0,0,0.1);
            }

            /* DataTable active row styling */
            table.dataTable tr.active td,
            table.dataTable tr.active {
                box-shadow: inset 0 0 0 9999px #D1DBE5 !important;
                color: black !important;
            }

            /* Download section */
            .download-wrapper { max-width: 1200px; margin: auto; padding: 20px; }
            .file-table {
                width: 100%; table-layout: fixed; border-collapse: collapse; margin-top: 20px;
            }
            .file-table th, .file-table td {
                width: 30%; padding: 10px; text-align: center;
                border-bottom: 1px solid #eaeaea;
            }
            .file-table th { background-color: #f5f5f5; font-weight: bold; }
            .download-button { padding: 6px 12px; font-size: 14px; }

            /* Disabled picker items */
            .bootstrap-select .disabled { opacity: 0.4 !important; }
            
            .file-table {
                table-layout: fixed;
                width: 100%;
            }
            .file-table th:nth-child(1),
            .file-table td:nth-child(1) {   /* File Name */
                width: 20%;
            }
            .file-table th:nth-child(2),
            .file-table td:nth-child(2) {   /* Description */
                width: 45%;
            }
            .file-table th:nth-child(3),
            .file-table td:nth-child(3) {   /* Size */
                width: 15%;
            }
            .file-table th:nth-child(4),
            .file-table td:nth-child(4) {   /* Action */
                width: 15%;
            }

            /* Footer styling */
            .footer {
                background-color: #193D5B; 
                color: white; 
                text-align: center;
                padding: 10px; 
            }
        "))
    ),
    
    navbarPage(
        title = tags$a(
            tagList(icon("project-diagram"), "The Single-Cell Hypertension Atlas"),
            href  = "#shiny-tab-home",
            style = "text-decoration: none; color: inherit; font-size: 22px;"
        ),
        id = "navbar",
        
        # --- Home Tab ---
        tabPanel(
            title = tagList(icon("house"), "Home"), value = "shiny-tab-home",
            
            # Intro header
            tags$head(tags$style(HTML("
                .header-intro {
                background-image: url('background.png');
                background-size: cover;
                background-position: center;
                text-align: center;
                padding: 100px 40px;
                margin: -20px -15px 20px -15px;
                position: relative; overflow: hidden;
                }
                .header-intro::before {
                content: ''; position: absolute; inset: 0; z-index: 0;
                }
                .header-intro > h1, .header-intro > p {
                position: relative; z-index: 1; color: #193D5B;
                }
                .header-intro > h1 { font-size: 42px; font-weight: bold; }
                .header-intro > p { font-size: 18px; }
            "))),
            div(
                class = "header-intro",
                h1("The Single-Cell Hypertension Atlas"), br(),
                p("The Single-Cell Hypertension Atlas offers an interactive platform to explore and analyze single-cell data in hypertension, supporting discovery and data sharing within the research community.",
                  style = "font-size: 20px; color: #193D5B;")
            ),
            
            # Study Design section
            div(
                class = "container-main", style = "margin-top: 40px; border-radius: 4px;",
                tags$h4(icon("clipboard", class = "fas"), tags$strong("Study Design", style = "font-size: 22px;")),
                p(  
                    "This study used three animal models to investigate hypertension, collecting five key organs for single-nucleus RNA and multiome sequencing. ",
                    "The schematic shows strain-specific treatments, tissue sampling, and assay types with corresponding sample sizes.",
                    style = "font-size: 16px; color: #2b2b2b; margin-bottom: 20px;"
                ),
                div(
                    style = "text-align: center;", 
                    tags$img(
                        src = "mouse.png", 
                        style = "max-width: 85%;", 
                        alt = "Tissue Image"
                        )
                    )
            ),
            
            # Clustering section
            div(
                class = "container-main", style = "margin-top: 40px; border-radius: 4px;",
                tags$h4(icon("layer-group", class = "title-icon"), tags$strong("Clustering Across Strains and Tissues", style = "font-size: 22px;")),
                div(style = "text-align: center;", tags$img(src = "tissue.png", style = "max-width: 88%;", alt = "Tissue Image")),
                tags$p(tags$strong("Tissue Abbreviations:"), style = "font-size: 16px; margin-top: 10px;"),
                div(
                    style = "display: flex; justify-content: center;",
                    tags$table(style = "border-collapse: collapse; font-size: 15px;", tags$tbody(
                        tags$tr(
                            tags$td("• HYP (hypothalamus)", style = "padding: 4px 16px;"), 
                            tags$td("• LV (left ventricle)", style = "padding: 4px 16px;"), 
                            tags$td("• LK (left kidney)", style = "padding: 4px 16px;"),
                            tags$td("• MCA (middle cerebral artery)", style = "padding: 4px 16px;"),
                            tags$td("• MSA (3rd mesenteric artery)", style = "padding: 4px 16px;")
                        )
                    ))
                ),
                tags$p(tags$strong("Cell Type Abbreviations:"), style = "font-size: 16px; margin-top: 10px;"),
                div(
                    style = "display: flex; justify-content: center;",
                    tags$table(style = "border-collapse: collapse; font-size: 15px;", tags$tbody(
                        tags$tr(
                            tags$td("• OPC (oligodendrocyte progenitor cell)"), tags$td("• NFO (newly formed oligodendrocyte)"),
                            tags$td("• CM (cardiomyocyte)"), tags$td("• POD (podocyte)")
                        ),
                        tags$tr(
                            tags$td("• PT (proximal tubule)"), tags$td("• TL (thin limb)"),
                            tags$td("• TAL (thick ascending limb)"), tags$td("• DCT (distal convoluted tubule)")
                        ),
                        tags$tr(
                            tags$td("• CT (connecting tubule)"), tags$td("• CD (collecting duct)"),
                            tags$td("• IC (intercalated cell)"), tags$td("• EC (endothelial cell)")
                        ),
                        tags$tr(
                            tags$td("• VSMC (vascular smooth muscle cell)"), tags$td("• DC (dendritic cell)"),
                            tags$td("• NK (natural killer cell)"), tags$td("• NKT (natural killer T cell)")
                        ),
                    ))
                ),
                br()
            ),
            
            # Paper Citation
            div(
                class = "container-main", style = "margin-top: 40px; border-radius: 4px;",
                tags$h4(icon("bookmark", class = "fas"), tags$strong("Paper Citation", style = "font-size: 22px;")),
                p("If you wish to cite this web resource or the associated dataset, please consider citing:", style = "font-size: 16px;")
            ),
            br()
            
        ),
        
        # --- DEG Tab ---
        tabPanel(
            title = tagList(icon("bars"), "DEG"),
            fluidPage(
                div(
                    class = "container-main", style = "margin-top: 40px; border-radius: 4px;",
                    
                    # Filters
                    h3(icon("dna"), "Differentially Expressed Genes List", style = "color: #193D5B; font-weight: bold;"),
                    p(
                        "This table displays differentially expressed genes (DEGs) filtered by strain, tissue, comparison group, and cell type. ",
                        "You can sort or search the table, and select any row to view corresponding gene expression plots below.",
                        style = "font-size: 16px;"
                    ),
                    fluidRow(column(
                        12,
                        div(
                            style = "background-color: #f9f9f9; border: 1px solid #ddd; border-radius: 6px; padding: 15px; margin-bottom: 10px;",
                            fluidRow(
                                column(2, pickerInput("strain",     "Strain",    choices = NULL)),
                                column(2, pickerInput("tissue",     "Tissue",    choices = NULL)),
                                column(3, pickerInput("comparison", "Comparison",choices = NULL)),
                                column(
                                    3, pickerInput(
                                        "celltype", "Cell Type", choices = NULL, multiple = TRUE,
                                        options = list(
                                            `actions-box`          = TRUE,
                                            `selected-text-format` = "count > 3",
                                            `dropup-auto`          = FALSE,
                                            `close-on-select`      = FALSE
                                        )
                                    )
                                ),
                                column(2, div(style = "margin-top: 25px;",
                                              actionButton("reset_filters", "Reset Filters", icon("redo"), class = "btn btn-primary btn-block")))
                            )
                        )
                    )),
                    br(),
                    
                    # DEG table & download
                    DT::dataTableOutput("deg_table"),
                    fluidRow(column(12, align = "right", downloadButton("download_deg", "Download DEGs", class = "btn btn-primary"))),
                    tags$hr(style = "border-top: 2px solid #ccc; margin: 20px 0;"),
                    
                    # Gene Expression Plots
                    h3(icon("images", class = "far"), "Gene Expression", style = "color: #193D5B; font-weight: bold;"),
                    p(
                        "After selecting a gene from the table above, this section will display:",
                        tags$ul(
                            tags$li("A UMAP plot showing spatial expression of the selected gene across all cells."),
                            tags$li("Two violin plots showing expression levels across cell types under the control and treatment conditions.")
                        ),
                        style = "font-size: 16px;"
                    ),
                    fluidRow(
                        column(4, plotOutput("umap_plot")),
                        column(8, plotOutput("violin_plot"))
                    ),
                    br(),
                    br()
                )
            )
        ),
        
        # --- Gene Tab ---
        tabPanel(
            title = tagList(icon("image"), "Gene"),
            fluidPage(
                div(
                    class = "container-main", style = "margin-top: 40px; border-radius: 4px;",
                    
                    h3(icon("image"), "Gene Expression Visualization", style = "color: #193D5B; font-weight: bold;"),
                    p(
                        "Use the search box below to select a gene of interest. This module displays expression and differential expression profiles across different cell types, tissues, strains, and treatments.",
                        style = "font-size: 16px;"
                    ),
                    tags$ul(
                        style = "font-size: 15px; padding-left: 20px;",
                        tags$li("The left plot shows overall gene expression (log2(expression+1)) and the percentage of expressing cells."),
                        tags$li("The right plot shows significant differential expression results (log2 fold change, p-value) from Wilcoxon tests."),
                        tags$li("Scroll below to view and download the underlying expression data table for the selected gene.")
                    ),
                    
                    div(
                        style = "background-color: #f9f9f9; border: 1px solid #ddd; border-radius: 6px; padding: 15px; margin-bottom: 10px;",
                        pickerInput("gene_select", "Search your gene:", choices = NULL,
                                    options = pickerOptions(liveSearch = TRUE, size = 10))
                    ), 
                    br(),
                
                    div(
                        style = "width: 100%; overflow-x: auto; text-align: center;",
                        div(style = "display: inline-block; min-width: 800px;",
                            plotOutput("gene_plot", width = "100%", height = "800px"))
                    ), br(),
                    
                    tags$hr(style = "border-top: 2px solid #ccc; margin: 20px 0;"),
                    h4("Gene Expression Table", style = "color: #193D5B; font-weight: bold;"),
                    p("This table displays the gene expression data used to generate the dot plot above for the selected gene.", style = "font-size: 16px;"),
                    br(),
                    DT::dataTableOutput("gene_table"),
                    fluidRow(column(12, align = "right", downloadButton("download_expression_table", "Download Expression Table", class = "btn btn-primary"))),
                    br()
                )
            )
        ),
        
        # --- Download Tab ---
        tabPanel(
            title = tagList(icon("cloud-arrow-down"), "Download"),
            fluidPage(
                div(
                    class = "container-main", style = "margin-top: 40px; border-radius: 4px;",
                    h3(icon("download"), "Download Center", style = "color: #193D5B; font-weight: bold;"),
                    
                    # Seurat Objects
                    div(
                        class = "download-wrapper",
                        h3("Seurat Objects (12 files)"),
                        p("Download Seurat objects for all conditions and tissues."),
                        tags$table(
                            class = "file-table",
                            tags$thead(tags$tr(tags$th("File Name"), tags$th("Description"), tags$th("Size"), tags$th("Action"))),
                            tags$tbody(
                                lapply(seurat_files, function(fname) {
                                    tags$tr(
                                        tags$td(fname),
                                        tags$td(seurat_descriptions[[fname]] %||% ""), 
                                        tags$td(sprintf("%.2f MB", file.info(file.path("www/seurat_objects", fname))$size / 1024^2)),
                                        tags$td(downloadButton(paste0("dl_", fname), "Download"))
                                    )
                                })
                            )
                        ),
                        
                        br(),
                        br(),
                        
                        # DEG Lists
                        h3("DEG Lists (2 files)"),
                        p("Download differentially expressed gene lists."),
                        tags$table(
                            class = "file-table",
                            tags$thead(tags$tr(tags$th("File Name"), tags$th("Description"), tags$th("Size"), tags$th("Action"))),
                            tags$tbody(
                                lapply(deg_files, function(fname) {
                                    tags$tr(
                                        tags$td(fname),
                                        tags$td(deg_descriptions[[fname]] %||% ""), 
                                        tags$td(sprintf("%.2f MB", file.info(file.path("www/deg_download", fname))$size / 1024^2)),
                                        tags$td(downloadButton(paste0("dl_", fname), "Download"))
                                    )
                                })
                            )
                        )
                    )
                )
            )
        )
    ),
    
    # Footer
    div(class = "footer", 
        p("Citation: ",
        "Some reference or attribution here.")
        )
)

# End of UI definition


# server.R

library(shiny)
library(DT)
library(shinyWidgets)
library(Seurat)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(ggsci)
library(ggh4x)
library(patchwork)

# ---- 1. Working directory & Data loading ----
setwd("/xdisk/mliang1/jingh/Rshiny/Shiny_Hypertension")

# DEG table
deg_data <- read.csv("www/deg_list/DEG.list.filtered.csv")
colnames(deg_data) <- c(
    "gene_symbol","cell_type","tissue","p_val","log2fc(ctrl vs trt)","pct.1","pct.2","p_val_adj","strain","control","treatment","comparison" 
    )


# Strain → tissue/comparison map
strain_map <- list(
    "C57BL/6" = list(
        tissue     = c("HYP", "LV", "LK"),
        comparison = c("Saline 3d vs AngII 3d", "Saline 3d vs AngII 28d")
    ),
    "SD" = list(
        tissue     = c("HYP", "LV", "LK", "MSA"),
        comparison = c("LS vs HS 3d")
    ),
    "SS" = list(
        tissue     = c("HYP", "LV", "LK", "MSA"),
        comparison = c("LS vs HS 3d", "LS vs HS 21d")
    ),
    "WKY" = list(
        tissue     = c("HYP", "LK", "MSA"),
        comparison = c("10w vs 26w")
    ),
    "SHR" = list(
        tissue     = c("HYP", "MCA", "LK", "LV", "MSA"),
        comparison = c("10w vs 26w")
    )
)

# Build celltype_map[strain][tissue][comparison]
celltype_map <- list()
for (st in unique(deg_data$strain)) {
    celltype_map[[st]] <- list()
    df_st <- filter(deg_data, strain == st)
    for (ti in unique(df_st$tissue)) {
        df_ti <- filter(df_st, tissue == ti)
        celltype_map[[st]][[ti]] <- list()
        for (cmp in unique(df_ti$comparison)) {
            df_cmp <- filter(df_ti, comparison == cmp)
            celltype_map[[st]][[ti]][[cmp]] <- sort(unique(df_cmp$cell_type))
        }
    }
}

# All unique values
all_strains     <- sort(unique(deg_data$strain))
all_tissues     <- sort(unique(deg_data$tissue))
all_comparisons <- sort(unique(deg_data$comparison))
all_celltypes   <- sort(unique(deg_data$cell_type))


# ---- 2. Load Seurat objects for UMAP/Violin ----
seurat_objects <- list(
    "C57BL/6_HYP" = readRDS("www/seurat_object/C57BL6_HYP.rds"),
    "C57BL/6_LV"  = readRDS("www/seurat_object/C57BL6_LV.rds"),
    "C57BL/6_LK"  = readRDS("www/seurat_object/C57BL6_LK.rds"),
    "SP_HYP"      = readRDS("www/seurat_object/SP_HYP.rds"),
    "SP_LK"       = readRDS("www/seurat_object/SP_LK.rds"),
    "SP_LV"       = readRDS("www/seurat_object/SP_LV.rds"),
    "SP_MCA"      = readRDS("www/seurat_object/SP_MCA.rds"),
    "SP_MSA"      = readRDS("www/seurat_object/SP_MSA.rds"),
    "SS_HYP"      = readRDS("www/seurat_object/SS_HYP.rds"),
    "SS_LK"       = readRDS("www/seurat_object/SS_LK.rds"),
    "SS_LV"       = readRDS("www/seurat_object/SS_LV.rds"),
    "SS_MSA"      = readRDS("www/seurat_object/SS_MSA.rds")
)


# ---- 3. Data for “Gene” tab ----
source("00.initial_setting.R")

deg_merged <- read.table(
    "www/deg_list/deg_download/DEG.all.out",
    sep = "\t", header = TRUE
)

expr_all <- bind_rows(
    read.table(
        "www/deg_list/deg_download/pseudo.DEG.all.out",
        sep = "\t", header = TRUE
    ) %>% select(
        pct = pct.1, avg_expr = avg_expr.1,
        gene_name, cell_type, project, strain, tissue,
        treatment = control
    ),
    read.table(
        "www/deg_list/deg_download/pseudo.DEG.all.out",
        sep = "\t", header = TRUE
    ) %>% select(
        pct = pct.2, avg_expr = avg_expr.2,
        gene_name, cell_type, project, strain, tissue,
        treatment
    )
)

# Preserve original cell order from initial settings
deg_merged$cell_type <- factor(deg_merged$cell_type, levels = cell_order)
expr_all$cell_type  <- factor(expr_all$cell_type,  levels = cell_order)

# Global theme
theme_set(theme_classic(base_size = 12))


# ---- 4. Server logic ----
server <- function(input, output, session) {
    
    # 4.1 Reactive filters & reset flag
    filters <- reactiveValues(
        strain     = "All",
        tissue     = "All",
        comparison = "All",
        celltype   = "All"
    )
    is_resetting <- reactiveVal(FALSE)
    
    # 4.2 Sync picker inputs → filters
    observeEvent(input$strain,     { filters$strain     <- input$strain })
    observeEvent(input$tissue,     { filters$tissue     <- input$tissue })
    observeEvent(input$comparison, { filters$comparison <- input$comparison })
    
    # 4.3 “All” logic for multi-select celltype
    observeEvent(input$celltype, {
        isolate({
            sel <- input$celltype
            if (is.null(sel) || length(sel) == 0) sel <- "All"
            if ("All" %in% sel && length(sel) > 1) {
                if (tail(sel, 1) == "All") {
                    sel <- "All"
                } else {
                    sel <- setdiff(sel, "All")
                }
            }
            updatePickerInput(session, "celltype", selected = sel)
            filters$celltype <- sel
        })
    })
    
    # 4.4 Compute valid options based on filters
    valid_options <- reactive({
        df <- deg_data
        if (filters$strain     != "All") df <- df[df$strain     == filters$strain, ]
        if (filters$tissue     != "All") df <- df[df$tissue     == filters$tissue, ]
        if (filters$comparison != "All") df <- df[df$comparison == filters$comparison, ]
        if (!("All" %in% filters$celltype)) {
            df <- df[df$cell_type %in% filters$celltype, ]
        }
        list(
            strains     = sort(unique(df$strain)),
            tissues     = sort(unique(df$tissue)),
            comparisons = sort(unique(df$comparison)),
            celltypes   = sort(unique(df$cell_type))
        )
    })
    
    # 4.5 Update all pickers with gray-out logic
    observe({
        opts <- valid_options()
        
        # Strain
        updatePickerInput(
            session, "strain",
            choices    = c("All" = "All", setNames(all_strains, all_strains)),
            selected   = filters$strain,
            choicesOpt = list(
                disabled = c(
                    FALSE,
                    sapply(all_strains, function(s) {
                        sub <- deg_data
                        if (filters$tissue     != "All") sub <- sub[sub$tissue     == filters$tissue, ]
                        if (filters$comparison != "All") sub <- sub[sub$comparison == filters$comparison, ]
                        !(s %in% sub$strain)
                    })
                )
            )
        )
        
        # Tissue
        updatePickerInput(
            session, "tissue",
            choices    = c("All" = "All", setNames(all_tissues, all_tissues)),
            selected   = filters$tissue,
            choicesOpt = list(
                disabled = c(
                    FALSE,
                    sapply(all_tissues, function(ti) {
                        sub <- deg_data
                        if (filters$strain     != "All") sub <- sub[sub$strain     == filters$strain, ]
                        if (filters$comparison != "All") sub <- sub[sub$comparison == filters$comparison, ]
                        !(ti %in% sub$tissue)
                    })
                )
            )
        )
        
        # Comparison
        updatePickerInput(
            session, "comparison",
            choices    = c("All" = "All", setNames(all_comparisons, all_comparisons)),
            selected   = filters$comparison,
            choicesOpt = list(
                disabled = c(
                    FALSE,
                    sapply(all_comparisons, function(cmp) {
                        sub <- deg_data
                        if (filters$strain != "All") sub <- sub[sub$strain == filters$strain, ]
                        if (filters$tissue != "All") sub <- sub[sub$tissue == filters$tissue, ]
                        !(cmp %in% sub$comparison)
                    })
                )
            )
        )
        
        # Celltype
        valid_ct <- if (
            filters$strain != "All" &&
            filters$tissue != "All" &&
            filters$comparison != "All"
        ) {
            tryCatch(
                celltype_map[[filters$strain]][[filters$tissue]][[filters$comparison]],
                error = function(e) character(0)
            )
        } else {
            all_celltypes
        }
        
        sel_ct <- filters$celltype
        if (is.null(sel_ct) || length(sel_ct) == 0) sel_ct <- "All"
        enabled <- union(valid_ct, sel_ct)
        disabled_flags <- !(all_celltypes %in% enabled)
        
        updatePickerInput(
            session, "celltype",
            choices    = c("All" = "All", setNames(all_celltypes, all_celltypes)),
            selected   = sel_ct,
            choicesOpt = list(disabled = c(FALSE, disabled_flags))
        )
    })
    
    
    # ---- 5. DEG table & download ----
    
    # Reactive subset
    filtered_deg <- reactive({
        req(input$strain, input$tissue, input$comparison, input$celltype)
        df <- deg_data
        if (input$strain     != "All") df <- df[df$strain     == input$strain, ]
        if (input$tissue     != "All") df <- df[df$tissue     == input$tissue, ]
        if (input$comparison != "All") df <- df[df$comparison == input$comparison, ]
        if (!("All" %in% input$celltype)) {
            df <- df[df$cell_type %in% input$celltype, ]
        }
        df
    })
    
    deg_proxy <- dataTableProxy("deg_table")
    
    # Render table
    output$deg_table <- renderDataTable({
        df <- filtered_deg()
        if (nrow(df) == 0) return(NULL)
        
        df <- df %>%
            mutate(
                p_val       = format(p_val, scientific = TRUE, digits = 3),
                `log2fc(ctrl vs trt)`  = format(`log2fc(ctrl vs trt)`, scientific = TRUE, digits = 3),
                p_val_adj   = format(p_val_adj, scientific = TRUE, digits = 3)
            )
        
        datatable(
            df,
            style     = "bootstrap",
            selection = "single",
            options   = list(
                dom         = 'Blfrtip',
                scrollX     = TRUE,
                buttons     = c("csv"),
                filter      = 'top',
                pageLength  = 10,
                lengthMenu  = c(10, 20, 50),
                columnDefs  = list(
                    list(visible = FALSE, targets = which(colnames(df) == "comparison") - 1)
                )
            ),
            rownames = FALSE
        )
    })
    
    # Download DEGs
    output$download_deg <- downloadHandler(
        filename = function() {
            paste0("filtered_DEG_", Sys.Date(), ".csv")
        },
        content = function(file) {
            write.csv(filtered_deg(), file, row.names = FALSE)
        }
    )
    
    
    # ---- 6. Reset Filters ----
    observeEvent(input$reset_filters, {
        is_resetting(TRUE)
        updatePickerInput(session, "strain",     choices = c("All" = "All", setNames(all_strains, all_strains)), selected = "All")
        updatePickerInput(session, "tissue",     choices = c("All" = "All", setNames(all_tissues, all_tissues)), selected = "All")
        updatePickerInput(session, "comparison", choices = c("All" = "All", setNames(all_comparisons, all_comparisons)), selected = "All")
        updatePickerInput(session, "celltype",   choices = c("All" = "All", setNames(all_celltypes, all_celltypes)), selected = "All")
        
        filters$strain     <- "All"
        filters$tissue     <- "All"
        filters$comparison <- "All"
        filters$celltype   <- "All"
        
        selectPage(deg_proxy, 1)
        invalidateLater(1000, session)
        observe({ is_resetting(FALSE) })
    })
    
    
    # ---- 7. UMAP & Violin for selected gene ----
    
    # Map strain → prefix for object lookup
    strain_to_prefix <- list(
        "C57BL/6" = "C57BL/6",
        "SD"       = "SS",
        "SS"       = "SS",
        "SHR"      = "SP",
        "WKY"      = "SP"
    )
    
    # Which row is selected?
    selected_deg_row <- reactive({
        req(input$deg_table_rows_selected)
        filtered_deg()[input$deg_table_rows_selected, ]
    })
    
    # Choose correct Seurat object & subset
    selected_seurat <- reactive({
        deg <- selected_deg_row()
        key <- paste0(strain_to_prefix[[deg$strain]], "_", deg$tissue)
        subset(seurat_objects[[key]], subset = strain == deg$strain)
    })
    
    # UMAP
    output$umap_plot <- renderPlot({
        deg <- selected_deg_row()
        FeaturePlot(
            selected_seurat(),
            features  = deg$gene_symbol,
            reduction = "umap"
        ) + scale_colour_gradientn(
            colours = c("#F1FAEE", "#A8DADC", "#457B9D", "#1D3557") 
        )+
            ggtitle(deg$gene_symbol) + labs(x = "UMAP 1", y = "UMAP 2")
    })
    
    # Violin
    output$violin_plot <- renderPlot({
        deg <- selected_deg_row()
        seu <- selected_seurat()
        ctrl <- deg$control
        trt  <- deg$treatment
        
        seu_sub <- subset(seu, subset = treatment %in% c(ctrl, trt))
        VlnPlot(
            seu_sub,
            features  = deg$gene_symbol,
            group.by  = "subclass_level1",
            split.by  = "treatment",
            pt.size   = 0.1
        ) +
            ggtitle(paste("Expression of", deg$gene_symbol, "in", deg$tissue, " (", ctrl, "vs", trt, ")")) +
            scale_x_discrete(expand = expansion(mult = c(0.05, 0.05))) +
            scale_fill_manual(
                values = setNames(c("#457B9D", "#E76F51"), c(ctrl, trt))
            )+ labs(x= "")
    })
    
    
    # ---- 8. “Gene” tab: picker + plots + table + download ----
    
    # Initialize gene choices
    available_genes <- sort(unique(c(deg_merged$gene_name, expr_all$gene_name)))
    updatePickerInput(session, "gene_select", choices = available_genes, selected = "Npr1")
    
    # Combined plot
    output$gene_plot <- renderPlot({
        req(input$gene_select)
        g <- input$gene_select
        
        # DEG‐Wilcox
        df_w <- filter(deg_merged, !is.na(p_val), p_val < 0.05, gene_name == g)
        p_w <- if (nrow(df_w) > 0) {
            ggplot(df_w, aes(x = treatment, y = cell_type)) +
                geom_point(aes(fill = -avg_log2FC, size = -log10(p_val), shape = p_val_adj < 0.05)) +
                scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
                scale_shape_manual(values = c(23, 21)) +
                labs(
                    title = paste(g, "(DEG‐wilcox)"),
                    fill  = "log2(FC)", size = "-log10(p‐value)", shape = "adj.p<0.05"
                ) +
                facet_nested(tissue ~ project + strain, scales = "free", space = "free") +
                theme(axis.text.x = element_text(angle = 45, hjust = 1))
        }
        
        # Overall expression
        df_e <- filter(expr_all, gene_name == g)
        p_e <- ggplot(df_e, aes(x = treatment, y = cell_type)) +
            geom_point(aes(color = log2(avg_expr + 1), size = pct * 100)) +
            scale_color_gradient(low = "white", high = "red") +
            labs(
                title = paste(g, "(overall expression)"),
                color = "log2(expr+1)", size = "Pct cells"
            ) +
            facet_nested(tissue ~ project + strain, scales = "free", space = "free") +
            theme(axis.text.x = element_text(angle = 45, hjust = 1))
        
        # Combine
        if (!is.null(p_w)) {
            (p_e + p_w) + plot_layout(widths = c(1.5, 1))
        } else {
            p_e
        }
    }, height = 800, width = 1000)
    
    # Expression table
    output$gene_table <- renderDataTable({
        req(input$gene_select)
        datatable(
            filter(expr_all, gene_name == input$gene_select),
            style     = "bootstrap",
            selection = "single",
            options   = list(pageLength = 10, scrollX = TRUE)
        )
    })
    
    # Download expression table
    output$download_expression_table <- downloadHandler(
        filename = function() {
            paste0("expression_data_", input$gene_select, ".csv")
        },
        content = function(f) {
            write.csv(
                filter(expr_all, gene_name == input$gene_select),
                f, row.names = FALSE
            )
        }
    )
    
    
    # ---- 9. Download tab handlers ----
    
    # Seurat object downloads
    seurat_files <- list.files("www/seurat_object", full.names = FALSE)
    for (fn in seurat_files) {
        local({
            file <- fn
            output[[paste0("dl_", file)]] <- downloadHandler(
                filename = function() file,
                content  = function(dest) {
                    file.copy(file.path("www/seurat_object", file), dest)
                }
            )
        })
    }
    
    # DEG list downloads
    deg_files <- list.files("www/deg_list/deg_download", full.names = FALSE)
    for (fn in deg_files) {
        local({
            file <- fn
            output[[paste0("dl_", file)]] <- downloadHandler(
                filename = function() file,
                content  = function(dest) {
                    file.copy(file.path("www/deg_list/deg_download", file), dest)
                }
            )
        })
    }
    }
    


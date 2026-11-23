group_review_ui <- function(id) {
  ns <- shiny::NS(id)
  shiny::tagList(
    shiny::uiOutput(ns("status")),
    shiny::conditionalPanel(
      "output.ready", ns = ns,
      shiny::tabsetPanel(
        id = ns("view"),
        shiny::tabPanel(
          "Suggestions", value = "suggestions",
          shiny::tags$div(
            class = "review-controls",
            shiny::selectInput(
              ns("filter"), "Show lines",
              c("All flagged lines" = "all", "With a suggestion" = "suggested",
                "Needs review" = "unresolved"), width = "190px"
            ),
            shiny::numericInput(ns("page"), "Page", 1, min = 1, step = 1, width = "90px")
          ),
          shiny::uiOutput(ns("row_status")),
          shiny::tags$div(class = "haplotype-qc-table", shiny::tableOutput(ns("suggestions")))
        ),
        shiny::tabPanel(
          "UMAP", value = "umap",
          shiny::actionButton(ns("build_umap"), "Build UMAP", class = "btn-preview"),
          shiny::uiOutput(ns("map"))
        )
      ),
      shiny::tags$section(
        class = "review-inspector",
        shiny::selectizeInput(
          ns("line"), "Inspect a germplasm line", choices = NULL,
          options = list(placeholder = "Search by germplasmName"), width = "100%"
        ),
        shiny::uiOutput(ns("line_details")),
        shiny::uiOutput(ns("decision_controls")),
        shiny::uiOutput(ns("decision_download"))
      )
    )
  )
}

group_review_server <- function(id, current_report, settings_pending, report) {
  shiny::moduleServer(id, function(input, output, session) {
    ns <- session$ns
    page_size <- 20L
    review <- shiny::reactive(current_report()$review)
    map_result <- shiny::reactiveVal(NULL)
    batch_ids <- shiny::reactiveVal(character())
    decisions <- shiny::reactiveVal(NULL)

    select_line <- function(ids, selected) {
      shiny::updateSelectizeInput(session, "line", choices = c("", ids), selected = selected, server = TRUE)
    }

    shiny::observeEvent(report()$key, {
      map_result(NULL)
      batch_ids(character())
      decisions(NULL)
      result <- report()$review
      if (is.null(result)) return()
      ids <- rownames(result$raw_calls)
      focused <- result$suggestions$germplasmName[1L]
      if (is.na(focused)) focused <- ""
      select_line(ids, focused)
    }, ignoreNULL = FALSE)

    output$ready <- shiny::reactive(!is.null(review()))
    shiny::outputOptions(output, "ready", suspendWhenHidden = FALSE)

    output$status <- shiny::renderUI({
      if (settings_pending()) {
        return(shiny::helpText("Run Haplotype Analysis to review groups with the current settings."))
      }
      result <- review()
      if (is.null(result)) return(shiny::helpText("Run Haplotype Analysis to review group assignments."))
      suggestions <- result$suggestions
      suggested <- sum(!is.na(suggestions$Suggested_Group))
      shiny::tags$p(
        role = "status",
        sprintf("%d lines flagged at a split SNP. %d have a suggested group and %d need review.",
                nrow(suggestions), suggested, nrow(suggestions) - suggested),
        "See Group Review Suggestions in the Excel report."
      )
    })

    filtered_suggestions <- shiny::reactive({
      rows <- review()$suggestions
      if (is.null(rows)) return(NULL)
      if (length(batch_ids())) rows <- rows[rows$germplasmName %in% batch_ids(), , drop = FALSE]
      if (identical(input$filter, "suggested")) rows <- rows[!is.na(rows$Suggested_Group), , drop = FALSE]
      if (identical(input$filter, "unresolved")) rows <- rows[is.na(rows$Suggested_Group), , drop = FALSE]
      rows
    })

    shiny::observeEvent(list(input$filter, batch_ids(), review()), {
      pages <- max(1L, ceiling(NROW(filtered_suggestions()) / page_size))
      shiny::updateNumericInput(session, "page", value = 1, max = pages)
    })

    page_range <- shiny::reactive({
      total <- NROW(filtered_suggestions())
      if (!total) return(integer())
      page <- input$page
      if (is.null(page) || !is.finite(page)) page <- 1L
      page <- min(max(1L, floor(page)), ceiling(total / page_size))
      seq.int((page - 1L) * page_size + 1L, min(page * page_size, total))
    })

    output$row_status <- shiny::renderUI({
      rows <- page_range()
      shiny::tagList(
        if (length(batch_ids())) shiny::tags$p(
          sprintf("Showing the %d flagged lines selected on the map. ", length(batch_ids())),
          shiny::actionButton(ns("clear_batch"), "Clear selection", class = "btn-xs")
        ),
        shiny::helpText(if (length(rows)) {
          sprintf("Rows %d–%d of %d", rows[1L], utils::tail(rows, 1L), NROW(filtered_suggestions()))
        } else {
          "No flagged lines match this view."
        })
      )
    })

    output$suggestions <- shiny::renderTable({
      shiny::req(length(page_range()))
      rows <- filtered_suggestions()[page_range(), , drop = FALSE]
      data.frame(
        germplasmName = rows$germplasmName,
        "Original group" = rows$Original_Group,
        "Suggested group" = ifelse(is.na(rows$Suggested_Group), "Needs review", rows$Suggested_Group),
        "Runner-up" = rows$Runner_Up_Group,
        "Supporting SNPs" = rows$Supporting_Markers,
        Missing = rows$Missing,
        Heterozygous = rows$Explicit_Heterozygous,
        Unresolved = rows$Unresolved,
        check.names = FALSE
      )
    }, rownames = FALSE, na = "")

    selected_suggestion <- shiny::reactive({
      rows <- review()$suggestions
      index <- match(shiny::req(input$line), rows$germplasmName)
      if (!is.na(index)) rows[index, , drop = FALSE]
    })

    output$line_details <- shiny::renderUI({
      line <- shiny::req(input$line)
      result <- shiny::req(review())
      assignment <- current_report()$assignments
      index <- match(line, assignment$germplasmName)
      shiny::req(!is.na(index))
      detail <- result$evidence(line)
      suggestion <- selected_suggestion()
      typed <- detail$Call_Type %in% c("Homozygous", "Heterozygous")
      marker_order <- order(!detail$Distinct_Signatures, !typed, is.na(detail$Best_Allele))
      shown <- utils::head(detail[marker_order, , drop = FALSE], 30L)

      bases <- c("A", "C", "G", "T")
      allele_chip <- function(label) {
        base <- toupper(label)
        style <- if (base %in% bases) paste0(" review-allele-", base) else ""
        shiny::tags$span(class = paste0("review-allele", style), label)
      }
      allele <- function(value) {
        if (is.na(value) || !nzchar(value)) return(allele_chip("–"))
        parts <- trimws(strsplit(value, "[/|*]")[[1L]])
        if (length(parts) == 2L && all(toupper(parts) %in% bases)) {
          shiny::tagList(allele_chip(parts[1L]), "/", allele_chip(parts[2L]))
        } else {
          allele_chip(value)
        }
      }
      candidate_label <- function(column, default) {
        group <- shown[[column]][1L]
        if (is.null(group) || is.na(group)) default else paste("Group", group)
      }
      shiny::tagList(
        shiny::tags$p(
          shiny::tags$strong(line),
          paste(" • Original group", assignment$Concatenated_Grouping[index]),
          if (!is.null(suggestion)) {
            if (is.na(suggestion$Suggested_Group)) " • Needs review" else
              paste(" • Suggested group", suggestion$Suggested_Group)
          }
        ),
        if (!is.null(suggestion)) shiny::tagList(
          shiny::helpText(suggestion$Status),
          if (is.finite(suggestion$Best_Score)) shiny::helpText(
            sprintf("Weighted agreement on %d SNPs (%.0f%% of those that distinguish these groups). Group %s %.0f%%. Group %s %.0f%%.",
                    suggestion$Comparable_Markers, 100 * suggestion$Comparable_Coverage,
                    suggestion$Best_Group, 100 * suggestion$Best_Score,
                    suggestion$Runner_Up_Group, 100 * suggestion$Runner_Up_Score)
          ),
          shiny::helpText(
            sprintf("%d SNPs favor the best group and %d favor the runner-up. Among QC-passed markers, %d calls are missing, %d heterozygous, and %d unresolved.",
                    suggestion$Supporting_Markers, suggestion$Conflicting_Markers, suggestion$Missing,
                    suggestion$Explicit_Heterozygous, suggestion$Unresolved)
          )
        ),
        if (nrow(shown)) shiny::tags$div(
          class = "haplotype-qc-table",
          shiny::tags$table(
            class = "table table-condensed",
            shiny::tags$caption(sprintf("SNP letters and group signatures. Showing %d of %d markers.",
                                        nrow(shown), nrow(detail))),
            shiny::tags$thead(shiny::tags$tr(lapply(
              c("SNP", "Line call", candidate_label("Best_Group", "Best match"),
                "Evidence", candidate_label("Runner_Up_Group", "Runner-up"), "Evidence"),
              shiny::tags$th
            ))),
            shiny::tags$tbody(lapply(seq_len(nrow(shown)), function(i) {
              shiny::tags$tr(
                shiny::tags$td(shown$Marker[i]),
                shiny::tags$td(allele(shown$Raw_Call[i])),
                shiny::tags$td(allele(shown$Best_Allele[i])),
                shiny::tags$td(shown$Best_Result[i]),
                shiny::tags$td(allele(shown$Runner_Up_Allele[i])),
                shiny::tags$td(shown$Runner_Up_Result[i])
              )
            }))
          )
        ) else shiny::helpText("No group signature is available for this line.")
      )
    })

    output$decision_controls <- shiny::renderUI({
      row <- selected_suggestion()
      if (is.null(row)) return(NULL)
      saved <- decisions()
      saved_row <- match(row$germplasmName, saved$germplasmName)
      selected <- if (!is.na(saved_row)) saved$Reviewed_Group[saved_row] else row$Suggested_Group
      if (is.na(selected)) selected <- ""
      groups <- sort(unique(current_report()$assignments$Concatenated_Grouping), method = "radix")
      shiny::tags$div(
        class = "review-controls",
        shiny::selectInput(
          ns("decision_group"), "Your reviewed group",
          c("Keep for review" = "", stats::setNames(groups, groups)), selected = selected, width = "210px"
        ),
        shiny::actionButton(ns("save_decision"), "Save Review Decision", class = "btn-primary"),
        if (!is.na(saved_row)) shiny::tags$span("Decision saved for this run.")
      )
    })

    shiny::observeEvent(input$save_decision, {
      row <- shiny::req(selected_suggestion())
      shiny::req(!is.null(input$decision_group))
      group <- input$decision_group
      shiny::req(!nzchar(group) || group %in% current_report()$assignments$Concatenated_Grouping)
      decision <- data.frame(
        germplasmName = row$germplasmName, Original_Group = row$Original_Group,
        Suggested_Group = row$Suggested_Group, Reviewed_Group = group,
        Decision = if (!nzchar(group)) "Needs review" else if (group == row$Original_Group) {
          "Confirmed original"
        } else {
          "Reassigned"
        }
      )
      saved <- decisions()
      if (NROW(saved)) saved <- saved[saved$germplasmName != row$germplasmName, , drop = FALSE]
      decisions(rbind(saved, decision))
    })

    output$decision_download <- shiny::renderUI({
      if (!NROW(decisions())) return(NULL)
      shiny::downloadButton(ns("download_decisions"), "Download Review Decisions")
    })

    output$download_decisions <- shiny::downloadHandler(
      filename = function() "group_review_decisions.csv",
      content = function(file) {
        shiny::req(review(), NROW(decisions()))
        utils::write.csv(decisions(), file, row.names = FALSE, na = "")
      }
    )

    shiny::observeEvent(input$build_umap, {
      result <- shiny::req(review())
      assignments <- current_report()$assignments
      map_lines <- data.frame(
        germplasmName = assignments$germplasmName,
        Group = assignments$Concatenated_Grouping,
        Needs_Review = assignments$germplasmName %in% result$suggestions$germplasmName
      )
      min_comparable <- max(2L, ceiling(ncol(result$raw_calls) / 2))
      map <- shiny::withProgress(message = "Building UMAP", value = 0.3, {
        tryCatch(group_review_umap(result$raw_calls, map_lines, min_comparable = min_comparable), error = function(e) {
          list(error = conditionMessage(e))
        })
      })
      map_result(map)
    })

    map_colors <- shiny::reactive({
      points <- shiny::req(map_result()$points)
      groups <- sort(unique(points$Group), method = "radix")
      stats::setNames(grDevices::hcl.colors(length(groups), "Dark 3"), groups)
    })

    output$map <- shiny::renderUI({
      map <- map_result()
      if (is.null(map)) return(shiny::helpText("Build UMAP to explore lines with enough comparable SNPs."))
      if (!is.null(map$error)) return(shiny::tags$p(role = "alert", class = "text-danger", map$error))
      excluded <- ngettext(map$excluded, "line excluded", "lines excluded")
      shiny::tagList(
        shiny::tags$p(sprintf("%d lines plotted using %d markers. %d %s.",
                              nrow(map$points), map$markers, map$excluded, excluded)),
        shiny::helpText("Colors show original groups. Outlined points are flagged at a split SNP. Click a point to inspect it or drag a box around lines to review together."),
        shiny::uiOutput(ns("map_legend")),
        shiny::plotOutput(
          ns("umap_plot"), height = "395px",
          click = ns("map_click"), hover = shiny::hoverOpts(ns("map_hover")),
          brush = shiny::brushOpts(ns("map_brush"), resetOnNew = TRUE)
        ),
        shiny::uiOutput(ns("map_hover_label")),
        shiny::uiOutput(ns("map_batch"))
      )
    })

    output$umap_plot <- shiny::renderPlot({
      points <- shiny::req(map_result()$points)
      colors <- map_colors()[points$Group]
      graphics::par(mar = c(4.2, 4.2, 0.5, 0.6), las = 1, bty = "l")
      graphics::plot(points$UMAP1, points$UMAP2, type = "n", xlab = "UMAP 1", ylab = "UMAP 2")
      graphics::points(points$UMAP1, points$UMAP2, pch = 16, cex = 0.9, col = colors)
      flagged <- points$Needs_Review
      graphics::points(points$UMAP1[flagged], points$UMAP2[flagged], pch = 21, cex = 1.3,
                       col = "#231b32", bg = colors[flagged], lwd = 1.4)
    })

    output$map_legend <- shiny::renderUI({
      colors <- map_colors()
      if (length(colors) > 16L) return(shiny::helpText("Hover over a line to see its group."))
      shiny::tags$div(class = "review-map-legend", lapply(names(colors), function(group) {
        shiny::tags$span(
          shiny::tags$span(class = "review-map-swatch", style = paste0("background:", colors[[group]])),
          group
        )
      }))
    })

    nearby_line <- function(position) {
      shiny::nearPoints(shiny::req(map_result()$points), shiny::req(position),
                        xvar = "UMAP1", yvar = "UMAP2", threshold = 14, maxpoints = 1)
    }

    output$map_hover_label <- shiny::renderUI({
      line <- nearby_line(input$map_hover)
      if (nrow(line)) shiny::helpText(sprintf("%s • Group %s%s", line$germplasmName, line$Group,
                                              if (line$Needs_Review) " • Flagged" else ""))
    })

    shiny::observeEvent(input$map_click, {
      line <- nearby_line(input$map_click)
      if (nrow(line)) select_line(rownames(review()$raw_calls), line$germplasmName)
    })

    map_batch <- shiny::reactive({
      points <- map_result()$points
      if (is.null(input$map_brush) || is.null(points)) return(character())
      chosen <- shiny::brushedPoints(points, input$map_brush, xvar = "UMAP1", yvar = "UMAP2")
      chosen$germplasmName[chosen$Needs_Review]
    })

    output$map_batch <- shiny::renderUI({
      selected <- map_batch()
      if (length(selected)) shiny::tags$div(
        class = "review-controls",
        shiny::tags$span(sprintf("%d flagged lines in the selection.", length(selected))),
        shiny::actionButton(ns("review_batch"), "Review Selected Lines")
      )
    })

    shiny::observeEvent(input$review_batch, {
      ids <- shiny::req(map_batch())
      batch_ids(ids)
      shiny::updateSelectInput(session, "filter", selected = "all")
      select_line(rownames(review()$raw_calls), ids[1L])
      shiny::updateTabsetPanel(session, "view", selected = "suggestions")
    })
    shiny::observeEvent(input$clear_batch, batch_ids(character()))
  })
}

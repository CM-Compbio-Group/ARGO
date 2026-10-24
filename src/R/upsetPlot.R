
export_width  = 3   # Match this with your pdf() or ggsave() width
export_height = 3   # Match this with your pdf() or ggsave() height

# 2. Allocate relative percentage fractions (e.g., 30% for left bars, 50% for top bars)
left_bars_fraction = 0.30  
top_bars_fraction  = 0.50  

calculated_left_width = export_width * left_bars_fraction
calculated_top_height = export_height * top_bars_fraction



# 1. Convert list data and generate combination matrix
mat = list_to_matrix(res.sig)
m = make_comb_mat(mat)

# 2. Get the full descending order indices based on size
ordered_set_indices = order(set_size(m_subset), decreasing = FALSE)


full_order = order(comb_size(m), decreasing = TRUE)

# Set up your relative height proportions (fractions between 0 and 1)
top_bars_proportion = 0.5   # The top bar plot will take up 50% of the total layout height
set_bars_proportion = 0.35  # The left bar plot will take up 35% of the total layout width


# 3. Subset the matrix to ONLY keep the top 10 intersections
m_subset = m[full_order[1:10]]
size_top_nums  = 8   # Controls numbers above the top bars & top axis ticks
size_top_title = 8   # Controls the vertical "Intersection size" axis title
size_set_names = 8   # Controls the right set names ("tradeSeq", "ARGO", etc.)
size_set_size_title = 8

# 4. Create the UpSet layout object using the subsetted matrix
p = UpSet(
  m_subset,
  pt_size = unit(2, "mm"),
  lwd = 0.75,
  comb_col = "black",   
  bg_col = "gray90",    
  
  # Tell it to strictly order the subset columns from 1 to 10 (descending)
  comb_order = seq_len(ncol(m_subset)), 
  set_order = ordered_set_indices,
  
  # Move horizontal bars to the left side
  right_annotation = NULL, 
  left_annotation = upset_left_annotation(
    m_subset, 
    width = unit(calculated_left_width, "inches"), 
    bar_width=0.4,
    annotation_name_gp = gpar(fontsize = size_set_size_title, fontfamily = "Arial"), 
    axis_param = list(gp = gpar(fontsize = size_set_names, fontfamily = "Arial"))
  ),
  
  # Set names font configurations
  row_names_gp = gpar(fontsize = size_set_names, fontfamily = "Arial"),
  
  # Top intersection size bar plot formatting
  top_annotation = upset_top_annotation(
    m_subset,
    height = unit(calculated_top_height, "inches"), 
    bar_width = 0.4, 
    add_numbers = TRUE, 
    numbers_gp = gpar(fontsize = size_top_nums, fontfamily = "Arial"),
    annotation_name_side = "left",
    annotation_name_rot = 90,
    annotation_name_gp = gpar(fontsize = size_top_title, fontfamily = "Arial"),
    axis_param = list(
      gp = gpar(fontsize = size_top_nums, fontfamily = "Arial") 
    )
  )
) 

pdf('figs/upsetPlot_FDR_05.pdf', width = export_width, height = export_height)
draw(p)
dev.off()


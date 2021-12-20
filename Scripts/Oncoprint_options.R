#### Oncoprint Options ####
cnnl_color <- rgb(244,221,81, maxColorValue = 255)
int_palette <- colorRampPalette(c(rgb(244,221,81, maxColorValue = 255),"firebrick2"))
ColorSchema <- c(
  "navyblue",
  "dodgerblue3",
  cnnl_color,
  int_palette(15)[8])
names(ColorSchema) <- c("HomoDel","HemiDel","CNNL","DelGain")
col = c(ColorSchema, SNV = "green4", Gain = "firebrick1", Uncertain_Homo = rgb(128, 128, 191, maxColorValue = 255), Uncertain_Gain = rgb(255, 130, 130, maxColorValue = 255), #Uncertain_Homo = "grey20", Uncertain_Gain = "grey30", 
        GrayZoneImb = "Gray50", GrayZoneNoImb = "Gray80", Germline = "springgreen")
col = c(col, Deletion = "lightseagreen")

get_type_fun = function(x) strsplit(x, ";")[[1]]

alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = "white", col = NA))
  },
  HomoDel = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col["HomoDel"], col = NA))
  },
  HemiDel = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col["HemiDel"], col = NA))
  },
  CNNL = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col["CNNL"], col = NA))
  },
  DelGain = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col["DelGain"], col = NA))
  },
  Gain = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col["Gain"], col = NA))
  },
  Uncertain_Gain = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
              gp = gpar(fill = col["Uncertain_Gain"], col = NA))
  },
  Uncertain_Homo = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
              gp = gpar(fill = col["Uncertain_Homo"], col = NA))
  },
  GrayZoneImb = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
              gp = gpar(fill = col["GrayZoneImb"], col = NA))
  },
  GrayZoneNoImb = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
              gp = gpar(fill = col["GrayZoneNoImb"], col = NA))
  },
  Germline = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.66,
              gp = gpar(fill = col["Germline"], col = NA))
  },
  Deletion = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col["Deletion"], col = NA))
  },
  SNV = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, 
              gp = gpar(fill = col["SNV"], col = NA))
  }
)

heatmap_legend_param = list(title = "Aberrations", at = c("HomoDel", "HemiDel", "CNNL", "DelGain", "Gain", "Germline", "SNV", "Deletion", 
                                                          "Uncertain_Homo", "Uncertain_Gain", "GrayZoneImb", "GrayZoneNoImb"), 
                            labels = c("HomoDel", "HemiDel", "CNNL", "Unb. Gain", "Bal. Gain", "Germline", "Non-Syn SNV", "Deletion on chrX", 
                                       "Uncertain_Homo", "Uncertain_Gain", "Imb_WTGrayZoneCN", "NoImb_WTGrayZoneCN"))
heatmap_legend_param = list(title = "Aberrations", at = c("HomoDel", "Uncertain_Homo", "HemiDel", "CNNL", "DelGain", "Uncertain_Gain", "Gain", 
                                                          "Germline", "SNV", "Deletion", "GrayZoneImb", "GrayZoneNoImb"), 
                            labels = c("HomoDel", "Uncertain_Homo", "HemiDel", "CNNL", "Unb. Gain", "Uncertain_Gain", "Bal. Gain", 
                                       "Germline", "Non-Syn SNV", "Deletion on chrX", "Imb_WTGrayZoneCN", "NoImb_WTGrayZoneCN"))

#### Oncoprint Options low TC ####
cnnl_color <- rgb(244,221,81, maxColorValue = 255)
int_palette <- colorRampPalette(c(rgb(244,221,81, maxColorValue = 255),"firebrick2"))
ColorSchema <- c(
  "navyblue",
  "dodgerblue3",
  cnnl_color,
  int_palette(15)[8])
names(ColorSchema) <- c("Deletion","HemiDel","CNNL","DelGain")
col_lowTC = c(ColorSchema, SNV = "green4", Gain = "firebrick1",
        GrayZoneImb = "Gray50", GrayZoneNoImb = "Gray80", Germline = "springgreen")
#col_lowTC = c(col, Deletion = "lightseagreen")

alter_fun_lowTC = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = "white", col = NA))
  },
  HemiDel = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col_lowTC["HemiDel"], col = NA))
  },
  CNNL = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col_lowTC["CNNL"], col = NA))
  },
  DelGain = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col_lowTC["DelGain"], col = NA))
  },
  GrayZoneImb = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
              gp = gpar(fill = col_lowTC["GrayZoneImb"], col = NA))
  },
  GrayZoneNoImb = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
              gp = gpar(fill = col_lowTC["GrayZoneNoImb"], col = NA))
  },
  Germline = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.66,
              gp = gpar(fill = col_lowTC["Germline"], col = NA))
  },
  Gain = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col_lowTC["Gain"], col = NA))
  },
  Deletion = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col_lowTC["Deletion"], col = NA))
  },
  SNV = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, 
              gp = gpar(fill = col_lowTC["SNV"], col = NA))
  }
)

heatmap_legend_param_lowTC = list(title = "Aberrations", at = c("Deletion", "HemiDel", "CNNL", "DelGain", "Gain", "Germline", "SNV", 
                                                          "GrayZoneImb", "GrayZoneNoImb"), 
                            labels = c("Likely Loss", "HemiDel", "CNNL", "Unb. Gain", "Likely Gain", "Germline", "Non-Syn SNV",
                                       "Imb_WTGrayZoneCN", "NoImb_WTGrayZoneCN"))

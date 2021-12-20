#### Oncoprint Options ####
cnnl_color <- rgb(244,221,81, maxColorValue = 255)
int_palette <- colorRampPalette(c(rgb(244,221,81, maxColorValue = 255),"firebrick2"))
ColorSchema <- c(
  "navyblue",
  "dodgerblue3",
  cnnl_color,
  int_palette(15)[8])
names(ColorSchema) <- c("HomoDel","HemiDel","CNNL","Unb.Gain")
col = c(ColorSchema, SNV = "green4", Bal.Gain = "firebrick1", "Uncertain HomoDel" = rgb(128, 128, 191, maxColorValue = 255), "Uncertain Bal.Gain" = rgb(255, 130, 130, maxColorValue = 255), #Uncertain_Homo = "grey20", Uncertain_Gain = "grey30", 
        Imb_WTGrayZoneCN = "Gray50", NoImb_WTGrayZoneCN = "Gray80", Germline = "springgreen")
col = c(col, "Deletion on chrX" = "lightseagreen")

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
  Unb.Gain = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col["Unb.Gain"], col = NA))
  },
  Bal.Gain = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col["Bal.Gain"], col = NA))
  },
  "Uncertain Bal.Gain" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
              gp = gpar(fill = col["Uncertain Bal.Gain"], col = NA))
  },
  "Uncertain HomoDel" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
              gp = gpar(fill = col["Uncertain HomoDel"], col = NA))
  },
  Imb_WTGrayZoneCN = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
              gp = gpar(fill = col["Imb_WTGrayZoneCN"], col = NA))
  },
  NoImb_WTGrayZoneCN = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
              gp = gpar(fill = col["NoImb_WTGrayZoneCN"], col = NA))
  },
  Germline = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.66,
              gp = gpar(fill = col["Germline"], col = NA))
  },
  "Deletion on chrX" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col["Deletion on chrX"], col = NA))
  },
  SNV = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, 
              gp = gpar(fill = col["SNV"], col = NA))
  }
)

# heatmap_legend_param = list(title = "Aberrations", at = c("HomoDel", "HemiDel", "CNNL", "DelGain", "Gain", "Germline", "SNV", "Deletion", 
#                                                           "Uncertain_Homo", "Uncertain_Gain", "GrayZoneImb", "GrayZoneNoImb"), 
#                             labels = c("HomoDel", "HemiDel", "CNNL", "Unb. Gain", "Bal. Gain", "Germline", "Non-Syn SNV", "Deletion on chrX", 
#                                        "Uncertain_Homo", "Uncertain_Gain", "Imb_WTGrayZoneCN", "NoImb_WTGrayZoneCN"))
heatmap_legend_param = list(title = "Aberrations", at = c("HomoDel", "Uncertain HomoDel", "HemiDel", "CNNL", "Unb.Gain", "Uncertain Bal.Gain", "Bal.Gain", 
                                                          "Germline", "SNV", "Deletion on chrX", "Imb_WTGrayZoneCN", "NoImb_WTGrayZoneCN"), 
                            labels = c("HomoDel", "Uncertain HomoDel", "HemiDel", "CNNL", "Unb.Gain", "Uncertain Bal.Gain", "Bal.Gain", 
                                       "Germline", "Non-Syn SNV", "Deletion on chrX", "Imb_WTGrayZoneCN", "NoImb_WTGrayZoneCN"))

#### Oncoprint Options low TC ####
cnnl_color <- rgb(244,221,81, maxColorValue = 255)
int_palette <- colorRampPalette(c(rgb(244,221,81, maxColorValue = 255),"firebrick2"))
ColorSchema <- c(
  "navyblue",
  "dodgerblue3",
  cnnl_color,
  int_palette(15)[8])
names(ColorSchema) <- c("Likely Loss","HemiDel","CNNL","Unb.Gain")
col_lowTC = c(ColorSchema, SNV = "green4", "Likely Gain" = "firebrick1",
        Imb_WTGrayZoneCN = "Gray50", NoImb_WTGrayZoneCN = "Gray80", Germline = "springgreen")
col_lowTC = c(col_lowTC, "Deletion on chrX" = "lightseagreen")

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
  Unb.Gain = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col_lowTC["Unb.Gain"], col = NA))
  },
  Imb_WTGrayZoneCN = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
              gp = gpar(fill = col_lowTC["Imb_WTGrayZoneCN"], col = NA))
  },
  NoImb_WTGrayZoneCN = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"),
              gp = gpar(fill = col_lowTC["NoImb_WTGrayZoneCN"], col = NA))
  },
  Germline = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.66,
              gp = gpar(fill = col_lowTC["Germline"], col = NA))
  },
  "Likely Gain" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col_lowTC["Likely Gain"], col = NA))
  },
  "Likely Loss" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col_lowTC["Likely Loss"], col = NA))
  },
  "Deletion on chrX" = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), 
              gp = gpar(fill = col["Deletion on chrX"], col = NA))
  },
  SNV = function(x, y, w, h) {
    grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, 
              gp = gpar(fill = col_lowTC["SNV"], col = NA))
  }
)

heatmap_legend_param_lowTC = list(title = "Aberrations", at = c("Likely Loss", "HemiDel", "CNNL", "Unb.Gain", "Likely Gain", "Germline", "SNV",
                                                          "Imb_WTGrayZoneCN", "NoImb_WTGrayZoneCN"), 
                            labels = c("Likely Loss", "HemiDel", "CNNL", "Unb.Gain", "Likely Gain", "Germline", "Non-Syn SNV",
                                       "Imb_WTGrayZoneCN", "NoImb_WTGrayZoneCN"))

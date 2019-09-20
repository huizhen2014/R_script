library(AnnotationForge)
makeOrgPackage(go=fGO,
               version="0.1",
               maintainer = "carlos <carlos@google.com>",
               author="carlos <carlos@google.com>",
               outputDir = ".",
               tax_id="58729",
               genus="Taeniopygia",
               species = "guttata",
               goTable = "go",verbose=T)

install.packages("./org.Tguttata.eg.db",repos=NULL,type="source")


#! /usr/bin/R

#############################
##REMOVED FROM THE ANALYSIS##
#############################

# These species were named by one or both of ALA and the Banksia phylogeny but
# were thrown out for reasons specified below:

# Species on ALA but removed from the analysis because they had less than 10
# occurrence points:
dspec <-
    c("Banksia_aculeata", "Banksia_anatona", "Banksia_aurantia", 
      "Banksia_bella", "Banksia_catoglypta", "Banksia_comosa",
      "Banksia_concinna", 
      "Banksia_corvijuga", "Banksia_croajingolensis", "Banksia_epica", 
      "Banksia_epimicta", "Banksia_foliolata", "Banksia_fuscobractea", 
      "Banksia_hirta", "Banksia_idiogenes", "Banksia_insulanemorecincta", 
      "Banksia_ionthocarpa", "Banksia_lepidorhiza", "Banksia_montana", 
      "Banksia_oligantha", "Banksia_prionophylla", "Banksia_pseudoplumosa", 
      "Banksia_recurvistylis", "Banksia_rosserae", "Banksia_rufistylis", 
      "Banksia_tricuspis", "Banksia_trifontinalis", "Banksia_viscida", 
      "Banksia_wonganensis")

# Species in the phylogeny but not on ALA:
# where cd = comparative.data(banksiaphy, banksiadata)
cd$drop[[1]][!(cd$drop[[1]] %in% dspec)]
[1] "Banksia_goodii"       "Banksia_dolichostyla"
# dolichostyla is meant to be a subsp of sphaerocarpa

# Species on ALA but not in the phylogeny:
nophy <- c("Banksia_acanthopoda", "Banksia_acuminata", "Banksia_cypholoba",
           "Banksia_erythrocephala", "Banksia_horrida", "Banksia_sclerophylla", 
           "Banksia_seneciifolia", "Banksia_xylothemelia")

######################
##SPECIES CATEGORIES##
######################

# East coast species:
ecs <-
    c("Banksia_aemula", "Banksia_aquilonia", "Banksia_canei",
      "Banksia_conferta", "Banksia_croajingolensis", "Banksia_dentata",
      "Banksia_ericifolia", "Banksia_integrifolia", "Banksia_marginata",
      "Banksia_oblongifolia", "Banksia_ornata", "Banksia_paludosa",
      "Banksia_plagiocarpa", "Banksia_robur", "Banksia_saxicola",
      "Banksia_serrata", "Banksia_spinulosa"
      )

# West coast species:
#wcs <- banksia.data$species[!(banksia.data$species %in% ecs)]

Dryandras <-
    c("Banksia_concinna", "Banksia_serra", "Banksia_idiogenes",
      "Banksia_glaucifolia", 
      "Banksia_ionthocarpa", "Banksia_biterax", "Banksia_formosa", 
      "Banksia_bella", "Banksia_foliosissima", "Banksia_mucronulata", 
      "Banksia_nobilis_subsp_fragrans", "Banksia_stuposa",
      "Banksia_nobilis_subsp_nobilis", 
      "Banksia_plumosa", "Banksia_pseudoplumosa", "Banksia_viscida", 
      "Banksia_anatona", "Banksia_falcata", "Banksia_prolata",
      "Banksia_tenuis_var_reptans", 
      "Banksia_tenuis_var_tenuifolia", "Banksia_obovata", "Banksia_stenoprion", 
      "Banksia_arctotidis", "Banksia_nivea_subsp_fuliginosa",
      "Banksia_nivea_subsp_nivea", 
      "Banksia_tortifolia", "Banksia_brunnea", "Banksia_dallaneyi", 
      "Banksia_fasciculata", "Banksia_densa", "Banksia_platycarpa", 
      "Banksia_insulanemorecincta", "Banksia_rufistylis", "Banksia_pallida", 
      "Banksia_comosa", "Banksia_rufa", "Banksia_kippistiana",
      "Banksia_columnaris", 
      "Banksia_splendida_subsp_splendida", "Banksia_nana",
      "Banksia_splendida_subsp_macrocarpa", 
      "Banksia_serratuloides_subsp_perissa", "Banksia_borealis",
      "Banksia_prionophylla", 
      "Banksia_tridentata", "Banksia_carlinoides", "Banksia_fraseri", 
      "Banksia_fuscobractea", "Banksia_calophylla", "Banksia_porrecta", 
      "Banksia_undata", "Banksia_arborea", "Banksia_corvijuga",
      "Banksia_heliantha", 
      "Banksia_epimicta", "Banksia_mimica", "Banksia_vestita",
      "Banksia_cynaroides", 
      "Banksia_sessilis_var_cordata", "Banksia_hewardiana",
      "Banksia_sessilis_var_sessilis", 
      "Banksia_strictifolia", "Banksia_echinata", "Banksia_wonganensis", 
      "Banksia_squarrosa_subsp_argillacea",
      "Banksia_subpinnatifida_var_imberbis", 
      "Banksia_squarrosa_subsp_squarrosa", "Banksia_polycephala",
      "Banksia_subpinnatifida_var_subpinnatifida", 
      "Banksia_armata", "Banksia_cirsioides", "Banksia_shanklandiorum", 
      "Banksia_catoglypta", "Banksia_fililoba", "Banksia_drummondii", 
      "Banksia_alliacea", "Banksia_pellaeifolia", "Banksia_lepidorhiza", 
      "Banksia_pteridifolia", "Banksia_octotriginta", "Banksia_purdieana", 
      "Banksia_meganotia", "Banksia_proteiodes", "Banksia_shuttleworthiana", 
      "Banksia_obtusa", "Banksia_bipinnatifida")

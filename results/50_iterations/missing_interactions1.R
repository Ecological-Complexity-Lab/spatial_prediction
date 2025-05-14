# ---- missing interactions ----
# install.packages("rglobi")
library(rglobi)

# e.g. get all “eats” interactions where the source is "Panthera leo"
ints_euporbia <- get_interactions_by_taxa(
  sourcetaxon   = "euphorbia balsamifera"
)

# Take a look
view(ints_euporbia)

ints_camponotus <- get_interactions_by_taxa(
  sourcetaxon   = "camponotus"
)

# Take a look
view(ints_camponotus)

ints_feae <- get_interactions_by_taxa(
  sourcetaxon   = "feae"
)

# Take a look
view(ints_feae)

ints_camponotus <- ints_feae %>% filter(source_taxon_name == "Camponotus feae")
view(ints_camponotus)

ints_anaspis <- get_interactions_by_taxa(
  sourcetaxon   = "anaspis proteus"
)

# Take a look
view(ints_anaspis)

ints_urell <- get_interactions_by_taxa(
  sourcetaxon   = "urelliosoma"
)

# Take a look
view(ints_urell)

ints_euph <- get_interactions_by_taxa(
  sourcetaxon   = "euphorbia balsamifera"
)

# Take a look
view(ints_euph)

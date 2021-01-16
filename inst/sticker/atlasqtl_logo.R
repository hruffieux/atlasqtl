require(hexSticker)
require(dplyr)
require(ggpubr)


seed <- 123
set.seed(seed)

p <- 20
q <- 1000
df <- data.frame("snp" = 1:p, 
                 "hotspot_size" = q*rbeta(p, shape1 = 0.1, shape2 = 100))

logo <- ggdotchart(df, x = "snp", y = "hotspot_size",
                   color = ifelse(df$hotspot_size>1, "#E6AB02", "#66A61E"),
                sorting = "none",                       # Sort value in descending order
                add = "segments",                             # Add segments from y = 0 to dots
                rotate = FALSE,                                # Rotate vertically
                dot.size = 1,                                 # Large dot size
                ggtheme = theme_void()) +                       # ggplot2 theme
  theme_transparent() + theme(legend.position = "none") 
logo

dir.create("man/figures/", showWarnings = FALSE)

sticker(logo, 
        package="atlasqtl", 
        p_size=4.5, 
        s_x=0.975, 
        s_y=1, 
        s_width=1.7, 
        s_height=1.3,
        p_x = 0.6, 
        p_y = 1.4, 
        u_color = "white", 
        u_size = 1,
        h_fill="black", 
        h_color="grey",
        filename="man/figures/atlasqtl_logo.png",
        dpi = 1200)

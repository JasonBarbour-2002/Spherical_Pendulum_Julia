module My_Themes
    using GLMakie
    function my_dark_Theme()
        Theme(
            fontsize = 20,
            textcolor = :white, 
            Slider = (
                color_active="#589AFA",
                color_active_dimmed="#FF9300",
                color_inactive=:black 
            ), 
            Button = (
                buttoncolor=:black, 
                labelcolor_active=:black,
                labelcolor_hover=:white, 
                buttoncolor_active = "#FF9300", 
                buttoncolor_hover= "#EE8200"
            ), 
            Heatmap = (
                highclip = :black,
                colormap = cgrad(:hot, rev=true),
                alpha = 0.7
            )
        ) 
    end
    function custum_dark_Theme()
        T = merge(my_dark_Theme(), theme_dark())
        merge(T, theme_latexfonts())
    end

    function my_light_Theme()
        Theme(
            fontsize = 20, 
            Heatmap = (
                highclip = :white,
                colormap = :coolwarm,
                alpha = 0.9
            )
        )       
    end
    function custum_light_Theme()
        merge(theme_latexfonts(), my_light_Theme())
    end
end
##Kinetic Data Graphing Functions
##Author: Jenna Veenstra
##Date: July 5, 2024

#Packages required:
library(ggplot2)
library(SiZer)
library(ggforce)
library(gridExtra)
library(grid)
library(ggpmisc)
library(DataCombine)
library(rlist)
library(broom)
library(tidyverse)
library(SOfun)
library(xlsx)
library(pzfx)
library(plyr)
library(cowplot)
library(english)
library(utils)
library(boot)
library(nlraa)
library(minpack.lm)
library(EnvStats)

#Defining Michaelis-Menten function for substrate dependence experiments
M_M <- function(S, Vmax, Km) {
  Vmax * S / (Km + S)
}

#Determine the index of potential outliers for linear regression models
  #Uses Cook's Distance and DFBETAS to determine which points are potential outliers
get_outlier_index <- function(model){
  
  n = length(model$fitted.values)                               #Sample size
  
  #Determine the most influential point based on Cook's distance
  cooksD <- cooks.distance(model)                               #Calculates the Cook's distance (CD) for each point in the model
  order_cooksD <- order(cooksD, decreasing = TRUE)              #Orders the CDs so the most influential point is first
  threshold_cd <- 3 * mean(cooksD, na.rm = TRUE)                #Calculates the threshold for a CD to be considered influential
  potential_cook <- which(cooksD > threshold_cd)                #Determines which of the CDs are above the threshold
  outlier_cook <- intersect(order_cooksD, potential_cook)       #Uses the ordered CDs to order the potential CDs
  
  #Determine the most influential point based on DFBETAS method
  dfbetas <- abs(dfbetas(model))[,2]                            #Calculates the DFBETAS for the model, takes the absolute value, and grabs the values for the slope coefficient (assuming model is a linear regression)
  order_dfbetas <- order(dfbetas, decreasing = TRUE)            #Orders the DFBETAS so the most influential point is first
  threshold_betas <- 2/sqrt(n)                                  #Calculates the threshold for a DFBETAS to be considered influential
  potential_betas <- which(dfbetas > threshold_betas)           #Determines which of the DFBETAS are above the threshold
  outlier_betas <- intersect(order_dfbetas, potential_betas)    #Uses the ordered DFBETAS to order the potential DFBETAS
  
  #Grab the most influential point (CD is given priority)
  union <- union(outlier_cook, outlier_betas)                   #Combine the indices of the potential outliers
  if (any(union == 1) & model$model$y_mod[1] == 0){             #The first point of the fits are (0,0) unless there is background (y_mod != 0)
    union <- union[-which(union == 1)]                          #Remove index 1 if it's already the origin
  }
  outlier_index <- union[1]                                     #Define the outlier_index
  return(outlier_index)
}


#* Fit a piecewise linear model where the slope is 0 after the change point
#* The model is very simple: lm( y ~ pmin(x, cp)) where cp is the change point. This
#* gives standard regression up to the change point, and then no change after the change point, i.e. a slope of 0.
#* Of course, you need to iterate to find the actual change point.
#The comments* above and original code was accessed from http://www.stat.sfu.ca/~cschwarz/Stat-650/Notes/MyPrograms on 5/15/19 but the link no longer works
#Original comments will be marked with '*', my comments will not
piecewise.linear.simple0<- function (x, y, middle = 1) 
  #*  Fits a piecewise linear model where the slope is 0 after the change point
  #No modifications were needed for this function  
{
  piecewise.linear.likelihood <- function(alpha, x, y) {
    N_ss <<- length(x)
    fit <- lm(y ~ pmin(x,alpha) )
    Beta <- coefficients(fit)
    Mu <- Beta[1] + Beta[2] * pmin(x,alpha) 
    SSE <- sum(fit$residuals^2)
    sigma2 <- SSE/N_ss
    likelihood <- sum(dnorm(y, mean = Mu, sd = sqrt(sigma2), log=TRUE))
    return(likelihood)
  }
  r <- range(x)
  offset <- r * (1 - middle)/2
  low <- min(x) + offset
  high <- max(x) - offset
  temp <- optimize(piecewise.linear.likelihood, c(low, high), 
                   x = x, y = y, maximum = TRUE)
  return(temp$maximum)
}


piecewise.linear0 <- function (x, y, middle = 1, CI = FALSE, bootstrap.samples = 1000, 
                               sig.level = 0.05) {
  x_mod <- x
  y_mod <- y
  removed_x <- c()
  ticker <- 0
  origin_manual <<- FALSE
  repeat { #Added repeat for outlier detection
    alpha <- piecewise.linear.simple0(x_mod, y_mod, middle)
    model <- lm(y_mod ~ pmin(x_mod,alpha))
    r_squared <- summary(model)$r.squared
    if (!rm_outliers){break} #break out of the repeat if rm_outliers is FALSE
    
    #Use outlier detection function to get the index for the influential point
    outlier_index <<- get_outlier_index(model)
    if (length(removed_x) >= 3 | is.na(outlier_index) | r_squared > r_threshold){break} #Don't remove the influential point if more than three points have already been removed (for typical time courses of 5-6 pts)
                                                                                        #Don't remove the influential point if the fit R2 value is above the r_threshold
    if (length(outlier_index) >= 1){ #Indicates an influential point was found
      #If the outlier is the origin and it was included in the experiment, the problem could be background. 
      #So, manually changing the point to the origin then seeing if it still is the biggest outlier
      if (x_mod[outlier_index] == 0 & ticker == 0){
        y_mod[outlier_index] <- 0
        origin_manual <<- TRUE
        ticker <- ticker + 1
        
      } else {
        removed_x <- c(removed_x, x_mod[outlier_index]) #keep record of which x values were removed
        x_mod <- x_mod[-outlier_index]                  #Remove the x-value of the influential point
        y_mod <- y_mod[-outlier_index]                  #Remove the y-value of the influential point
        ticker <- ticker + 1
      }
    } else {
      break #If no influential points are found, break out of the repeat
    }
  }
  outlier_index <- c()
  for (i in removed_x){
    outlier_index <- c(outlier_index, which(x == i)) #create the outlier_index outside of the repeat for future use
  }
  #Create a list for the output of this function
  out <- NULL
  out$change.point <- alpha #breakpoint of the piecewise linear regression
  out$model <- model #the linear regression model
  out$x <- seq(min(x), max(x), length = 600)  #x values to be used for predictions (for graphing the curve)
  out$y <- predict(out$model, data.frame(x_mod = out$x), interval = "confidence", level = 0.95) #predicted values based on the model
  out$CI <- CI #Whether the confidence interval was included
  out$outlier_index <- outlier_index #The indices for points identified as influential
  class(out) <- "PiecewiseLinear"
  if (CI == FALSE) {return(out)} else { #calculate the confidence interval for breakpoint and slope
    #browser()
    data <- data.frame(x = x_mod, y = y_mod) #dataset of the x and y values after outliers are removed
    my.cp <- function(data, index) { #Function to be iterated through for bootstrapping
      x <- data[index, 1]
      y <- data[index, 2]
      cp <- piecewise.linear.simple0(x, y)
      model <- lm(y ~ pmin(x,cp) )
      out <- c(cp, model$coefficients[2] )
      return(out)
    }
    boot.result <- boot(data, my.cp, R = bootstrap.samples) #perform bootsrapping to find the confidence intervals
    out$intervals <- apply(boot.result$t, 2, quantile, probs = c(sig.level/2, 
                                                                 1 - sig.level/2)) #extract the intervals
    colnames(out$intervals) <- c("Change.Point", "Initial.Slope")
    out$CI <- t(out$CI)
    return(out) #return the out list
  }
}

y_by_function <- function(y_lim, round_digit){ #A function to find the best y-axis limits and breaks
  while (y_lim/round_digit < 4){
    round_digit <- round_digit/2
  }
  return(round_digit)
}

#arguments for Graphing_Kinetic_Experiments
  #update if you add arguments
arg_list <- c("my_data", "x", "y", "adj_y", "group_by", "points_to_remove", "time_units_factor", "graph_title", "y_title", "x_title", "num_reps", "y_lim", "y_by", "x_lim", "x_by", "base_file_name", "base_folder_name", "V_over_E_file_name", "substrate", "substrate_name", "avg_groupings", "exp", "replace_files", "data_frame_name", "filter", "z_threshold", "r_threshold", "rm_outliers")

#New Functions
#Assign arguments as global variables so they will continue through functions in functions
globalize_arguments <- function(argument_ch, argument_ob){
  #browser()
  assign(argument_ch, argument_ob, envir = globalenv())
}

#Removes from the global environment so it doesn't interfere after the pipeline is finished
remove_globalized_arguments <- function(argument_ob){
  #browser()
  rm(argument_ob)
}

#calculate the slope for group
calculate_slopes <- function(group, my_data, exp = exp){
  #Create variable with same naming system as the group argument
  my_data$group <- paste0("rep_", my_data$replicate, "_", substrate_name, "_", my_data[,substrate])
  
  #Determine which indices the group is
  index <- which(my_data$group == group)
  #Grab the x and y variables for linear regression below
  x_var <- my_data[index,x]
  y_var <- my_data[index,adj_y]
  #Initialize the x and y variables that will have the influential points removed
  x_mod <- x_var
  y_mod <- y_var
  removed_x <- c()
  ticker <- 0
  origin_manual <<- FALSE
  if (sum(y_var, na.rm = TRUE) == 0){ #If all the points have already been filtered out (except (0,0)), then the group isn't used further
    slope_x <- NA
  } else {
    repeat{ #Same process as the repeat in piecewise.linear0()
      model <- lm(y_mod~x_mod)
      r_squared <- summary(model)$r.squared
      if (!rm_outliers){break}
      #outlier detection by Cook's Distance and DFBETAS
      outlier_index <<- get_outlier_index(model)
      if (length(x_mod) < 3){
        slope_x <- 0
        break
      }
      if (is.na(outlier_index) | r_squared > r_threshold){break}
      #If the outlier is the origin and it was included in the experiment, the problem could be background. 
      #So, manually changing the point to the origin then seeing if it still is the biggest outlier
      if (x_mod[outlier_index] == 0 & ticker == 0){
        y_mod[outlier_index] <- 0
        origin_manual <<- TRUE
        ticker <- ticker + 1
        break
      } else {
        removed_x <- c(removed_x, x_mod[outlier_index])
        x_mod <- x_mod[-outlier_index]
        y_mod <- y_mod[-outlier_index]
        ticker <- ticker + 1
      }
      for (i in removed_x) {
        rm_index <- index[which(my_data[index,x] == i)] #determine what indices in main dataset were removed above
        my_data$removed_factor[rm_index] <- factor(T, levels = c(F, T)) #Set the removed_factor for removed points to true for graphing later
        my_data[rm_index, adj_y] <- NA #set the adjusted y value to NA for removed points so it's not included in future fits/calculations
      }
    }
    slope_x <- summary(model)$coefficients[2,1]/time_units_factor #calculate the slope in terms of seconds
  }
  my_data$slope[index[1]] <- slope_x #Add slope for this group to main dataset
  
  ##To avoid many if else statements, always include DNA_Conc and Enz_Conc columns in my_data
  V_over_E <- (slope_x*my_data$DNA_Conc[index[1]])/my_data$Enz_Conc[index[1]] #Calculate the V/E value for M-M fits
  my_data$V_over_E[index[1]] <- V_over_E #Add V/E for this group to main dataset
  return(my_data) #return the main dataset so changes are 'saved'
}

#Graph one rep of my_data
graph_one_rep <- function(rep, my_data, num_reps)
  #rep is the replicate being graphed
  #my_data is the main dataset
  #num_reps is the max number of reps in my_data
  {
  #browser()
  if (file.exists(paste0(base_folder_name, base_file_name, my_data$protein[1], "_rep_", rep, ".png")) & replace_files == FALSE){
    print(paste0("rep_", rep, " already exists")) #If the file the graph will be written to already exists and replace_files isn't true, the function prints and ends
  } else {
    #Grab the data for replicate 'rep'
    graphing_data <- my_data[which(my_data$replicate == rep),]
    
    if (sum(!is.na(graphing_data[,y])) != 0){#if there is data for this rep
      #initialize the graph with the x and y variables
      graph <<- ggplot(graphing_data, aes(x = graphing_data[,x], y = graphing_data[,y])) 
      
      #When x_lim isn't given, find the max x value
      if (is.na(x_lim)){
        x_lim2 <- max(graphing_data$time_points, na.rm = T)
      } else {x_lim2 <- x_lim}
      #When y_lim isn't given, find the max y value and, using the y_by_function, find the best limit and breaks for the y axis
      if (is.na(y_lim)){
        y_lim2 <- max(graphing_data[,y], na.rm = T)
        num_zeros <- ifelse(y_lim2 >= 0.1, 0, str_count(gsub(".*[.]","",gsub("(?:(0+))[1-9].*","\\1",as.character(y_lim2))),"0"))
        round_digit <- as.numeric(paste0("0.", paste0(rep(0,num_zeros), collapse = ""), 1))
        round_digit <- y_by_function(y_lim2, round_digit)
        y_lim2 <- round_any(y_lim2, round_digit, f = ceiling)
        y_by2 <- round_digit  
      } else {
        y_lim2 <- y_lim
        y_by2 <- y_by
      }
      
      #Set up variables so each facet of the graph can be labeled correctly
      sub_vals <- unique(graphing_data[,substrate], na.rm = T)
      sub_pts <- graphing_data[,substrate]
        
      facet_labels <- c()
      
      #Some replicates vary enzyme concentration, others don't
      #Below adds enzyme concentration if it varies across the substrate concentrations
      for (i in sub_vals){
        enzyme_conc <- graphing_data$Enz_Conc[which(graphing_data[,substrate] == i)][1]
        if (length(unique(graphing_data$Enz_Conc)) == 1){
          facet_labels <- c(facet_labels, paste0(substrate_name, ": ", i))
        } else {
          facet_labels <- c(facet_labels, paste0(substrate_name, ": ", i, ", Enz Conc: ", enzyme_conc, " (nM)"))
        }
      }
      names(facet_labels) <- as.character(sub_vals)
        
      graphing_data$facet <- FALSE
      names(facet_labels) <- unique(paste0(sub_pts,graphing_data$facet))

      graph <<- graph + 
        #Have the graph include the time courses for each substrate concentration with the x-axis free, but y-axis is defined based on above calculations
        facet_wrap(.~fct_inorder(paste0(sub_pts,graphing_data$facet)), scales = "free", nrow = 2, labeller = as_labeller(facet_labels)) + scale_y_continuous(breaks = seq(from = 0, to = y_lim2, by = y_by2), limits = c(0, (y_lim2 + (y_lim2/16))), name = y_title) +
        
        #Add the data points with their shape based on the removed factor where points removed from the fit are 'X'
        geom_point(size = 2, aes(x = graphing_data[,x], y = graphing_data[,y], shape = removed_factor), stroke = as.numeric(graphing_data$removed_factor)) + scale_shape_manual(values = c(16,4), guide = "none") + 
        
        #Add linear regression lines with confidence intervals using the adjusted y values where influential points are removed. Label the graph with the R2 value
        stat_smooth(method = "lm", se = TRUE, fullrange = TRUE, aes(y = graphing_data[,adj_y]), na.rm = TRUE, alpha = 0.7) + stat_poly_eq(formula = y~x, parse = TRUE, aes(y = graphing_data[,adj_y], size = 3, fontface = "bold.italic")) + 
        
        #Add all the labels to the graph
        xlab(x_title) + ggtitle(label = graph_title) + theme(plot.title = element_text(hjust = 0.5, size = 36), axis.title = element_text(size = 30), axis.text = element_text(size = 15), strip.text = element_text(size = 15))
      
      if (!dir.exists(base_folder_name)){dir.create(base_folder_name, recursive = T)} #create the required directories if they don't already exist
      #save the graph using the the provided naming
      ggsave(plot = graph, filename = paste0(base_folder_name, base_file_name, my_data$protein[1], "_rep_", rep, ".png"), width = 15, height = 7, unit = "in", dpi = 300)
    } else { #if there isn't any data in graphing_data
      print(paste("Rep", rep, "is all NAs"))
    }
  }
}

graph_one_V_over_E <- function(rep, my_data, exp = exp)
  #rep is the replicate
  #my_data is the main dataset
  #exp is the experiment name
  ##See graph_one_rep for overlapping comments
  {
  #browser()
  if (file.exists(paste0(base_folder_name, base_file_name, my_data$protein[1], "_V_over_E_rep_", rep, ".png")) & replace_files == FALSE ){
    print(paste0("V_over_E_rep_", rep, " already exists"))
  } else {
    graphing_data <- my_data[which(my_data$replicate == rep),]
    
    #Change the provided graph_title to include V/E
    parts <- strsplit(graph_title, split = ": ")
    graph_title <- paste0(parts[[1]][1], ": V/E (s-1) (", parts[[1]][2], ")")
    
    #Grab the V/E values for the replicate (first index for each timecourse per calculate_slopes())
    fit_indices <- which(!is.na(graphing_data$V_over_E))
    fit_data <- graphing_data[fit_indices,]
    
    if (grepl("sub_dep",exp) & dim(fit_data)[1] > 1){ #When the experiment is a substrate dependence
      #set up initial guesses for the Michaelis-Menten fit and grab the V (y) and S (x) variables
      initial_guess <- list(Vmax = 0.5, Km = 50)
      V <- fit_data$V_over_E
      S <- fit_data[,substrate]
      
      #Initialize values for outlier detection
      x_mod <- S
      y_mod <- V
      removed_x <- c()
      ticker <- 0
      if (sum(y_mod, na.rm = TRUE) == 0){ #if there are no points then the kinetic parameters are NA
        kcat <- NA
        Km <- NA
        Rsquared <- NA
      } else {
        repeat{ #repeated until no more outliers are detected or the r_threshold is achieved
          #Use non-linear regression with the Levenberg-Marquardt adjustment
          model <- nlsLM(y_mod ~ M_M(x_mod, Vmax, Km), start = initial_guess)
          #find the residuals for the model
          residuals <- abs(resid(model))

          if (!rm_outliers){break} #Leave repeat if rm_outliers is FALSE
          #outlier detection using the MAD statistic
          mad_resid <- mad(residuals)
          order_mad <- order(residuals, decreasing = TRUE)
          threshold_mad <- 3*mad_resid
          potential_mad <- which(residuals > threshold_mad)
          outlier_index <- intersect(order_mad, potential_mad)[1]
          
          if (is.na(outlier_index) | Rsquared > r_threshold){break}
       
          removed_x <- c(removed_x, x_mod[outlier_index])
          x_mod <- x_mod[-outlier_index]
          y_mod <- y_mod[-outlier_index]
          ticker <- ticker + 1
        }
        #Modify graphing_data to indicate which values shouldn't be in the fit and be graphed as an 'X'
        for (i in removed_x) {
          rm_index <- fit_indices[which(graphing_data[fit_indices,substrate] == i)]
          graphing_data$removed_factor_V[rm_index] <- factor(T, levels = c(F, T))
          graphing_data$adj_V_over_E[rm_index] <- NA
          graphing_data$shape_V[rm_index] <- 4 
        }
        #Grab the kinetic parameters from the model
        kcat <- as.numeric(round(coef(model)["Vmax"],3))
        Km <- as.numeric(round(coef(model)["Km"],1))
        #Determine the confidence intervals for each kinetic paramter
        conf_int <- confint(model,level = 0.95)
        #Create an x-value series to predict y-values for graphing the curve
        newdata <- data.frame(x_mod = seq(0, max(S), length.out = 200))
        V_new <- predict(model, newdata)
        
        #Use bootstrapping to find the upper and lower limits of the models confidence interval
        df <- data.frame(x_mod = x_mod, y_mod = y_mod)
        boot_fn <- function(data, indices){
          sample_data <- data[indices, ]
          model <- nlsLM(y_mod ~ M_M(x_mod, Vmax, Km), data = sample_data, control = list(maxiter = 1000))
          new_pred <- predict(model, newdata = newdata)
          return(new_pred)}
        boot_results <- try(boot(data = df, statistic = boot_fn, R = 1000))
        
        if (inherits(boot_results, "try-error")){ #Possible for nlsLM() (less likely than nls()) to not converge during bootstrapping so 'catch' the error
          plot_data <- data.frame(x_mod = newdata, V_new = V_new)
          CI = FALSE
        } else {
          boot_ci <- t(sapply(1:nrow(newdata), function(i) {
            quantile(boot_results$t[,i], c(0.025,0.975))}))
          #If confidence intervals are determined, create plot_data for plotting the curve
          plot_data <- data.frame(x_mod = newdata, V_new = V_new, lwr = boot_ci[,1], upr = boot_ci[,2])
          CI = TRUE
        }
        #Add kcat and Km values to the graphing_data
        graphing_data$kcat <- kcat
        graphing_data$Km <- Km
      }
      
      #Create vectors for the table embedded on the graph
      kcat_table <- c("kcat (s-1)", "Fit" = kcat, "Lower CI" = round(conf_int[1,1],3), "Upper CI" = round(conf_int[1,2],3))
      Km_table <- c(paste0("Km", sub_units), "Fit" = Km, "Lower CI" = round(conf_int[2,1],3), "Upper CI" = round(conf_int[2,2],3))
    } else { #if experiment isn't substrate dependence, the table includes experimental conditions
      kcat_table <- c("DNA (nM)", "Enz (nM)", "ATP (uM)", "Mg (mM)")
      Km_table <- c(graphing_data$DNA_Conc[1],graphing_data$Enz_Conc[1],graphing_data$ATP_Conc[1], graphing_data$Mg_Conc[1], graphing_data$date_done[1])
    }
    #browser()
    if (sum(!is.na(graphing_data[,y])) != 0){
      sub_pts <- graphing_data[,substrate]
      
      y_lim2 <- y_lim
      y_by2 <- y_by
      x_lim2 <- x_lim
      
      if (is.na(x_lim)){x_lim2 <- max(graphing_data[,substrate], na.rm = T)}
      if (is.na(y_lim)){
        if(dim(plot_data)[2] == 2){
          y_lim2 <- max(graphing_data$V_over_E, kcat, na.rm = T)
        } else{y_lim2 <- max(plot_data$upr, graphing_data$V_over_E, kcat, na.rm = T)}
        num_zeros <- ifelse(y_lim2 >= 0.1, 0, str_count(gsub(".*[.]","",gsub("(?:(0+))[1-9].*","\\1",as.character(y_lim2))),"0"))
        round_digit <- as.numeric(paste0("0.", paste0(rep(0,num_zeros), collapse = ""), 1))
        round_digit <- y_by_function(y_lim2, round_digit)
        y_lim2 <- round_any(y_lim2, round_digit, f = ceiling)
        y_by2 <- round_digit
      }
      
      graph <- ggplot(graphing_data, aes(x = sub_pts, y = V_over_E)) + 
        #Add the all the y values but shape according to removed_factor_V
        geom_point(size = 2, shape = graphing_data$shape_V, stroke = as.numeric(graphing_data$removed_factor_V)) +
        #Setting of labels and axes
        labs(title = graph_title) + xlab(substrate_name) + ylab("V/E (s-1)") + theme(plot.title = element_text(hjust = 0.5, size = 12), axis.title = element_text(size = 11), axis.text = element_text(size = 10)) + ylim(c(0,y_lim2)) + xlim(c(0,x_lim2))

      if (grepl("sub_dep", exp) & dim(fit_data)[1] > 1){ #Use plot_data to graph the Michaelis-Menten curve
        if (CI == T) { #include the confidence interval in the plot
          graph <- graph + geom_line(data = plot_data, aes(x = x_mod, y = V_new), color = "purple") + geom_ribbon(data = plot_data, aes(x = x_mod, y = V_new, ymin = lwr, ymax = upr), fill = "grey", alpha = 0.7)
        } else { #Bootstrap failed so no confidence interval and red to indicate that
          graph <- graph + geom_line(data = plot_data, aes(x = x_mod, y = V_new), color = "red") + annotate("text", x = x_lim2/4, y = y_lim2/16, label = "CI couldn't be calculated", color = "red", size = 3, fontface = "bold")
        }
      }
      #Create the table to embed in the graph
      condition_table <- rbind(kcat_table, Km_table)
      rownames(condition_table) <- NULL
      
      #Add the table to the graph
      graph <- graph + annotation_custom(tableGrob(condition_table, theme = ttheme_default(base_size = 10, base_colour = "black")), xmin=(x_lim2/2), xmax=x_lim2, ymin=(0+y_lim2/30), ymax=y_lim2/2.5)
      

      if (!dir.exists(base_folder_name)){dir.create(base_folder_name, recursive = T)}
      ggsave(plot = graph, filename = paste0(base_folder_name, base_file_name, my_data$protein[1], "_V_over_E_rep_", rep, ".png"), width = 6, height = 4, unit = "in", dpi = 300)
    } else {
      print(paste("Rep", rep, "is all NAs"))
    }
    #Edit main dataset with all the changes made the graphing_data subset
    my_data[which(my_data$replicate == rep),] <- graphing_data
  }
  #Return main dataset so changes can be 'saved'
  return(my_data)
}


calculate_averages <- function(avg_across, graphing_data, exp = exp)
  #avg_across is the name for the group needing to be averaged
  #graphing_data is the inputted dataset
  #exp is the type of experiment
  {
  #browser()
  if (grepl("AST", exp)){ #Specific to active site titrations
    index <- which(graphing_data$group == avg_across)
    to_avg <- graphing_data[index,adj_y]
    data_index <- c()
    repeat {
      if (!rm_outliers){break}
      
      #Using the z-score to determine if a replicate value is an outlier
      fp_zscore <- scale(to_avg)
      if (!any(is.na(fp_zscore))){
        #Determine which zscores are above the threshold provided and if there are any, which is the max
        fp_out_index <- c()
        if(any(abs(fp_zscore) > z_threshold)){
          fp_out <- to_avg[which(abs(fp_zscore) > z_threshold)]
          out_index <- which(max(abs(fp_zscore)) == max(abs(fp_zscore)))
          #remove the point from to_avg
          to_avg <- to_avg[-out_index]
        } else {break}
      } else {break}
    }
    if (exists("fp_out")) { #keeping record of which values are removed from the average
      for (fp in fp_out) {
        fp_out_index <- c(fp_out_index, which(to_avg == fp))
        data_index <- c(data_index, which(graphing_data[,adj_y] == fp))
      }
    }
    
    #Adding record of above to graphing_data
    graphing_data[data_index,adj_y] <- NA
    graphing_data$shape_for_avg[data_index] <- 4
    graphing_data$removed_factor_for_avg[data_index] <- factor(T, levels = c(F, T))
    
    #calculate the average and standard deviation of the remaining to_avg values
    graphing_data$avg_prod_frac[index[1]] <- mean(to_avg, na.rm = T)
    graphing_data$sd_prod_frac[index[1]] <- ifelse(is.na(sd(to_avg, na.rm = T)), 0, sd(to_avg, na.rm = T))
    
  } else { #For non AST experiments (same overall structure)
    index <- which(graphing_data$group == avg_across)
    sub_conc <- strsplit(avg_across, "_")[[1]][2]
    to_avg_V_over_E <- graphing_data$V_over_E[index]
    data_index <- c()
    repeat{
      if (!rm_outliers){break}
      VE_zscore <- scale(to_avg_V_over_E)
      if (!any(is.na(VE_zscore))){
        VE_out_index <- c()
        if (any(abs(VE_zscore) > z_threshold)){
          VE_out <- to_avg_V_over_E[which(abs(VE_zscore) > z_threshold)]
          out_index <- which(max(abs(VE_zscore)) == max(abs(VE_zscore)))
          to_avg_V_over_E <- to_avg_V_over_E[-out_index]
        } else {break}
      } else {break}
    }
    
    if (exists("VE_out")) {
      for (VE in VE_out) {
        VE_out_index <- c(VE_out_index, which(to_avg_V_over_E == VE))
        data_index <- c(data_index, which(graphing_data$V_over_E == VE))
      }
    }
    graphing_data$adj_V_over_E_for_avg[data_index] <- NA
    graphing_data$shape_V_for_avg[data_index] <- 4
    graphing_data$removed_factor_V_for_avg[data_index] <- factor(T, levels = c(F, T))
    
    graphing_data$avg_V_over_E[index[1]] <- mean(to_avg_V_over_E, na.rm = T)
    graphing_data$sd_V_over_E[index[1]] <- ifelse(is.na(sd(to_avg_V_over_E, na.rm = T)), 0, sd(to_avg_V_over_E, na.rm = T))
  }
  #return the updated graphing data
  return(graphing_data)
}

graph_averages_V_over_E <- function(rep, graphing_data, exp = exp)
  #rep is the replicates included for the average
  #graphing_data is the inputted dataset
  #exp is the experiment name
  ##See graph_V_over_E for overlapping comments
  {
  
  #Collapse all the reps included into one character string
  reps <- paste0(rep, collapse = "_")
  if (file.exists(paste0(base_folder_name, V_over_E_file_name, graphing_data$protein[1], "_V_over_E_reps_avg_", reps, ".png", collapse = "_and_")) & replace_files == FALSE ){
    print(paste0("_V_over_E_reps_avg", reps, " already exists", collapse = "_and_"))
  } else {
    parts <- strsplit(graph_title, split = ": ")
    graph_title <- paste0(parts[[1]][1], ": V/E (s-1) (", parts[[1]][2], ")")
    fit_indices <- which(!is.na(graphing_data$avg_V_over_E))
    fit_data <- graphing_data[fit_indices,]

    if (grepl("sub_dep",exp) & dim(fit_data)[1] > 1){
      initial_guess <- list(Vmax = 0.5, Km = 50)
      V <- fit_data$avg_V_over_E
      S <- fit_data[,substrate]
      
      x_mod <- S
      y_mod <- V
      removed_x <- c()
      ticker <- 0
      if (sum(y_mod, na.rm = TRUE) == 0){
        kcat <- NA
        Km <- NA
        Rsquared <- NA
      } else {
        repeat{
          model <- nlsLM(y_mod ~ M_M(x_mod, Vmax, Km), start = initial_guess)
          residuals <- abs(resid(model))
          if (!rm_outliers) {break}
          mad_resid <- mad(residuals)
          order_mad <- order(residuals, decreasing = TRUE)
          threshold_mad <- 3*mad_resid
          potential_mad <- which(residuals > threshold_mad)
          outlier_index <- intersect(order_mad, potential_mad)[1]
          if (is.na(outlier_index) | Rsquared > r_threshold){break}
          removed_x <- c(removed_x, x_mod[outlier_index])
          x_mod <- x_mod[-outlier_index]
          y_mod <- y_mod[-outlier_index]
          ticker <- ticker + 1
        }
        
        for (i in removed_x) {
          rm_index <- fit_indices[which(graphing_data[fit_indices,substrate] == i)]
          graphing_data$removed_factor_avg[rm_index] <- factor(T, levels = c(F, T))
          graphing_data$avg_adj_V_over_E[rm_index] <- NA
          graphing_data$shape_avg[rm_index] <- 4 
        }

        kcat <- as.numeric(round(coef(model)["Vmax"],3))
        Km <- as.numeric(round(coef(model)["Km"],1))
        conf_int <- confint(model,level = 0.95)
        newdata <- data.frame(x_mod = seq(0, max(S), length.out = 200))
        V_new <- predict(model, newdata)
        df <- data.frame(x_mod = x_mod, y_mod = y_mod)
        boot_fn <- function(data, indices){
          sample_data <- data[indices, ]
          model <- nlsLM(y_mod ~ M_M(x_mod, Vmax, Km), data = sample_data, control = list(maxiter = 1000))
          new_pred <- predict(model, newdata = newdata)
          return(new_pred)}
        boot_results <- try(boot(data = df, statistic = boot_fn, R = 1000))
        
        if (inherits(boot_results, "try-error")){
          plot_data <- data.frame(x_mod = newdata, V_new = V_new)
          CI = FALSE
        } else {
          boot_ci <- t(sapply(1:nrow(newdata), function(i) {
            quantile(boot_results$t[,i], c(0.025,0.975))}))
          plot_data <- data.frame(x_mod = newdata, V_new = V_new, lwr = boot_ci[,1], upr = boot_ci[,2])
          CI = TRUE
        }
        
        graphing_data$kcat_avg <- kcat
        graphing_data$Km_avg <- Km
      }
      
      kcat_table <- c("kcat (s-1)", "Fit" = kcat, "Lower CI" = round(conf_int[1,1],3), "Upper CI" = round(conf_int[1,2],3))
      Km_table <- c(paste0("Km", sub_units), "Fit" = Km, "Lower CI" = round(conf_int[2,1],3), "Upper CI" = round(conf_int[2,2],3))
    } else {
      kcat_table <- c("DNA (nM)", "Enz (nM)", "ATP (uM)", "Mg (mM)")
      Km_table <- c(graphing_data$DNA_Conc[1],graphing_data$Enz_Conc[1],graphing_data$ATP_Conc[1], graphing_data$Mg_Conc[1], graphing_data$date_done[1])
    }
    #Ensuring that the shape is 'X' (4) for the removed_factor
    graphing_data$shape[which(graphing_data$removed_factor == TRUE)] <- 4
    shape_breaks <- unique(graphing_data$shape)
    if (any(shape_breaks == 4)){
      shape_breaks <- as.numeric(moveMe(shape_breaks, "4 first"))
      shape_labels <- c("Removed", unique(graphing_data$replicate))
    } else {shape_labels <- unique(graphing_data$replicate)}
    
    if (sum(!is.na(graphing_data[,y])) != 0){
      sub_pts <- graphing_data[,substrate]
      
      y_lim2 <- y_lim
      y_by2 <- y_by
      x_lim2 <- x_lim
      
      if (is.na(x_lim)){
        x_lim2 <- max(graphing_data[,substrate], na.rm = T)
      }
      if (is.na(y_lim)){
        if(dim(plot_data)[2] == 2){
          y_lim2 <- max(graphing_data$V_over_E, kcat, na.rm = T)
        } else{y_lim2 <- max(plot_data$upr, graphing_data$V_over_E, kcat, na.rm = T)}
        num_zeros <- ifelse(y_lim2 >= 0.1, 0, str_count(gsub(".*[.]","",gsub("(?:(0+))[1-9].*","\\1",as.character(y_lim2))),"0"))
        round_digit <- as.numeric(paste0("0.", paste0(rep(0,num_zeros), collapse = ""), 1))
        round_digit <- y_by_function(y_lim2, round_digit)
        y_lim2 <- round_any(y_lim2, round_digit, f = ceiling)
        y_by <- round_digit
      }
      graph <<- ggplot(graphing_data, aes(x = sub_pts)) + 
        #Titles and axes
        labs(title = graph_title) + scale_x_continuous(limits = c(0, x_lim2), name = substrate_name) + scale_y_continuous(breaks = seq(from = 0, to = y_lim2, by = y_by), limits = c(0, y_lim2), name = "V/E (s-1)") + theme(plot.title = element_text(hjust = 0.5, size = 12), axis.title = element_text(size = 11), axis.text = element_text(size = 10)) + 
        
        #The average values shaped according to whether they are included in the Michaelis-Menten fit
        geom_point(size = 3, aes(y = avg_V_over_E), shape = graphing_data$shape_avg, stroke = as.numeric(graphing_data$removed_factor_avg)) + 
        
        #Each replicates value colored by replicate and shaped according to whether they are included in calculating the average
        geom_point(aes(y = V_over_E, color = as.factor(replicate)), shape = graphing_data$shape_V_for_avg, stroke = as.numeric(graphing_data$removed_factor_V_for_avg), size = 2) + labs(color = "Replicate") + 
        
        #Add the errorbar based on standard deviation to the average points
        geom_errorbar(aes(ymin = avg_V_over_E - sd_V_over_E, ymax = avg_V_over_E + sd_V_over_E), width = 0.9, position = position_dodge(0.3))

      if (grepl("sub_dep",exp) & dim(fit_data)[1] > 1){
        if (CI == T) {
          graph <- graph + geom_line(data = plot_data, aes(x = x_mod, y = V_new), color = "purple") + geom_ribbon(data = plot_data, aes(x = x_mod, y = V_new, ymin = lwr, ymax = upr), fill = "grey", alpha = 0.7)
        } else {
          graph <- graph + geom_line(data = plot_data, aes(x = x_mod, y = V_new), color = "red") + annotate("text", x = x_lim2/4, y = y_lim2/16, label = "CI couldn't be calculated", color = "red", size = 3, fontface = "bold")
        }
      } 
      condition_table <- rbind(kcat_table, Km_table)
      rownames(condition_table) <- NULL
      
      graph <- graph + annotation_custom(tableGrob(condition_table, theme = ttheme_default(base_size = 9, base_colour = "black")), xmin=(x_lim2/4), xmax=x_lim2, ymin=0, ymax=y_lim2/5)

      
      if (!dir.exists(base_folder_name)){dir.create(base_folder_name, recursive = T)}
        ggsave(plot = graph, filename = paste0(base_folder_name, V_over_E_file_name, my_data$protein[1], "_V_over_E_reps_avg_", reps, ".png", collapse = "_and_"), width = 6, height = 4, unit = "in", dpi = 300)
      } else {
        print(paste("Rep", reps, "are all NAs"))
    }
  }
  return(graphing_data)
}

#Graph one replicate of Active Site Titration
graphing_one_rep_ast <- function(rep, my_data)
  #rep is the replicate being graphed
  #my_data is the main dataset
  {
  #Grab the data for replicate 'rep'
  graphing_data <- my_data[which(my_data$replicate == rep),]

  #Determine the index for all adjusted y values that aren't NA (already removed for fit)
  index <- which(!is.na(graphing_data[,adj_y]))
  
  #Fit the data to a piecewise linear regression using the piecewise.linear0() function (includes the Cook's distance and DFBETAS tests)
  sfit <- piecewise.linear0(graphing_data[index,x], graphing_data[index,adj_y], middle = 1, CI = T, bootstrap.samples = 1000, sig.level = 0.05)
  #create a dataset for plotting the regression
  fit_lines <- data.frame(x = unlist(sfit$x), y = unlist(sfit$y))
  #Extract the breakpoint and its confidence interval
  breakpoint <- round(sfit$change.point, 2)
  bp.down <- round(sfit$intervals[1,1], 1)
  bp.up <- round(sfit$intervals[2,1], 1)
  #Extract the Endpoint and R2 value
  endpoint <- round(max(fit_lines$y.fit, na.rm = T)*100, 1)
  Rsquared <- round(summary(sfit$model)$r.squared,3)
  outlier_index <- sfit$outlier_index
  graphing_data$removed_factor[outlier_index] <- as.factor(T)
    
  if (origin_manual == TRUE){#Set during the piecewise.linear0() function for if a non-zero value was an outlier at x = 0
    graphing_data[,y][which(graphing_data[,x] == 0)] <- 0
  }
  #Set the shape variable according to if a point was removed
  graphing_data$shape[which(graphing_data$removed_factor == TRUE)] <- 4
  graphing_data$shape[which(graphing_data$removed_factor == FALSE)] <- 16
  
  #vector to include in the return of this function
  calculated_values <- c(breakpoint, endpoint, Rsquared)

  
  if (file.exists(paste0(base_folder_name, base_file_name, my_data$protein[1], "_rep_", rep, ".png", collapse = "_and_")) & replace_files == FALSE ){
    print(paste0("_rep_", reps, " already exists", collapse = "_and_"))
    #Return the outliers and fit values for the regression fit
    return(list(outlier_index, calculated_values))
  } else { 
    #browser()
    graph <<- ggplot() + 
      #Add a line according to the regression fit and a confidence interval based on the bootstrapping in piecewise.linear0()
      geom_line(aes(x, y.fit), data = fit_lines, color = "black") + geom_ribbon(data = fit_lines, aes(x = x, ymin = y.lwr, ymax = y.upr), fill = "grey", alpha = 0.7) +
      
      #Add dashed lines for where the breakpoint and endpoint are
      geom_vline(xintercept = breakpoint, linetype = "dashed") + geom_hline(yintercept = endpoint/100, linetype = "dashed") + 
      
      #Add the data points shaped according to whether they were included in the fit
      geom_point(data = graphing_data, aes(x = graphing_data[,x], y = graphing_data[,y]), size = 2, shape = graphing_data$shape, stroke = as.numeric(graphing_data$removed_factor))
  }
  #Set the x and y axis limits if not already set
  if (is.na(x_lim)){
    x_lim <- max(graphing_data[,x])
    x_by <- x_lim/10
  }
  if (y_lim < max(fit_lines$y.upr)){
    y_lim <- max(fit_lines$y.upr)
  }
  graph <<- graph + 
    #Add the limits and breaks to the x and y axes and their titles
    scale_y_continuous(breaks = seq(from = 0, to = y_lim, by = y_by), limits = c(0,y_lim), name = y_title) + scale_x_continuous(breaks = seq(from = 0, to = x_lim, by = x_by), limits = c(0, x_lim), name = x_title) + 
    
    #Add the other labels
    ggtitle(label = graph_title) + theme(plot.title = element_text(hjust = 0.5, size = 12), axis.title = element_text(size = 11), axis.text = element_text(size = 10))
  
  #Create the table to be embedded in the graph
  labels_new <<- c("DNA (nM)", "Breakpoint", "CI", "Endpoint", "R^2")
  concentrations_new <<- c(graphing_data$DNA_Conc[1], breakpoint, paste0(bp.down, " - ", bp.up), paste(endpoint, "%"), Rsquared)
  condition_table <- cbind(labels_new, concentrations_new)
  colnames(condition_table) <- NULL
  
  graph <- graph + annotation_custom(tableGrob(condition_table, theme = ttheme_default(base_size = 10, base_colour = "black")), xmin=(x_lim/3)*2, xmax=x_lim, ymin=(0+y_lim/20), ymax=y_lim/3.5)
  if (!dir.exists(base_folder_name)){dir.create(base_folder_name, recursive = T)}
  ggsave(plot = graph, filename = paste0(base_folder_name, base_file_name, my_data$protein[1], "_rep_", rep, ".png"), width = 5, height = 4, unit = "in", dpi = 300)
  
  return(list(outlier_index, calculated_values))
}

graphing_averages_ast <- function(rep, graphing_data, exp = exp, grouping = grouping)
  #rep is the replicates included for the average
  #graphing_data is the inputted dataset
  #exp is the experiment name
  #grouping is the index for the set of replicates being averaged
  ##See graph_one_rep_ast for overlapping comments
  {
  
  #Set up variable names and adjusted values for removal
  avg_y <- "avg_prod_frac"
  avg_adj_y <- paste0("avg_", adj_y)
  graphing_data[,avg_adj_y] <- graphing_data[,avg_y]
  graphing_data[list_of_points_to_remove[[grouping]],avg_adj_y] <- NA
  reps <- paste0(rep, collapse = "_")

  #Grab the data to be included in the fit
  keep_points <- which(!is.na(graphing_data[,avg_adj_y]))
  fit_data <- graphing_data[keep_points,]

  #See graph_one_rep_ast
  index <- which(!is.na(fit_data[,avg_adj_y]))
  sfit <- piecewise.linear0(fit_data[index,x], fit_data[index,avg_adj_y], middle = 1, CI = T, bootstrap.samples = 1000, sig.level = 0.05)
  fit_lines <- data.frame(x = sfit$x, y = sfit$y)
  breakpoint <- round(sfit$change.point, 2)
  bp.down <- round(sfit$intervals[1,1], 1)
  bp.up <- round(sfit$intervals[2,1], 1)
  endpoint <- round(max(fit_lines$y.fit, na.rm = T)*100, 1)
  Rsquared <- round(summary(sfit$model)$r.squared,3)
  outlier_index <- sfit$outlier_index
  graphing_data$removed_factor_avg <- factor(F, levels = c(F,T))
  graphing_data$removed_factor_avg[outlier_index] <- as.factor(T)
  
  graphing_data$shape <- 16
  graphing_data$shape[which(graphing_data$removed_factor == TRUE)] <- 4
  if (file.exists(paste0(base_folder_name, V_over_E_file_name, graphing_data$protein[1], "_reps_avg_", reps, ".png")) & replace_files == FALSE ){
    print(paste0("_reps_avg_", reps, " already exists"))
  } else {  
    graphing_data$shape_avg <- 16
    graphing_data$shape_avg[which(graphing_data$removed_factor_avg == TRUE)] <- 4
        
    graph <<- ggplot() + 
      #Line and confidence interval for fit
      geom_line(data = fit_lines, aes(x, y.fit)) + geom_ribbon(data = fit_lines, aes(x = x, ymin = y.lwr, ymax = y.upr), fill = "grey", alpha = 0.7) + 
      
      #All the average values shaped based on what was included for the fit
      geom_point(data = graphing_data, size = 3, aes(x = graphing_data[,x], y = graphing_data[,avg_y]), shape = graphing_data$shape_avg, stroke = as.numeric(graphing_data$removed_factor_avg)) + 
      
      #All the individual replicate values colored by replicate and shaped based on whether they were included in the average
      geom_point(data = graphing_data, size = 2, aes(x = graphing_data[,x], y = graphing_data[,y], color = as.factor(replicate)), shape = graphing_data$shape, stroke = as.numeric(graphing_data$removed_factor)) +
      
      #Adding the error bar to the average values
      geom_errorbar(data = graphing_data, aes(x = graphing_data[,x], y = graphing_data[,avg_y], ymin = (graphing_data[,avg_y] - sd_prod_frac), ymax = (graphing_data[,avg_y] + sd_prod_frac)), width = max(graphing_data[,x]/20, na.rm = T)) + 
      
      #Adding the dashed lines indicating breakpoint and endpoint locations
      geom_vline(xintercept = breakpoint, linetype = "dashed") + geom_hline(yintercept = endpoint/100, linetype = "dashed")
      
    
    if (is.na(x_lim)){
      x_lim <- max(graphing_data[,x])
      if (is.na(x_by)){
        x_by <- x_lim/10
      }
    }
    if (y_lim < max(fit_lines$y.upr)){
      y_lim <- max(fit_lines$y.upr)
    }
    #Adding labels and axes limits
    graph <<- graph + scale_y_continuous(breaks = seq(from = 0, to = y_lim, by = y_by), limits = c(0,y_lim), name = y_title) + scale_x_continuous(breaks = seq(from = 0, to = x_lim, by = x_by), limits = c(0, x_lim), name = x_title) + 
      ggtitle(label = graph_title) + theme(plot.title = element_text(hjust = 0.5, size = 12), axis.title = element_text(size = 11), axis.text = element_text(size = 10), legend.position = "none") + scale_fill_discrete(guide = "none") #+labs(color = "Replicate")
    
    labels_new <<- c("DNA (nM)", "Breakpoint", "CI", "Endpoint", "R^2")
    concentrations_new <<- c(graphing_data$DNA_Conc[1], breakpoint, paste0(bp.down, " - ", bp.up), paste(endpoint, "%"), Rsquared)
    
    condition_table <- cbind(labels_new, concentrations_new)
    colnames(condition_table) <- NULL
    
    graph <- graph + annotation_custom(tableGrob(condition_table, theme = ttheme_default(base_size = 10, base_colour = "black")), xmin=(x_lim/3)*2, xmax=x_lim, ymin=(0+y_lim/20), ymax=y_lim/3.5)
    
    #browser()
    if (!dir.exists(base_folder_name)){dir.create(base_folder_name, recursive = T)}
    ggsave(plot = graph, filename = paste0(base_folder_name, V_over_E_file_name, graphing_data$protein[1], "_reps_avg_", reps, ".png"), width = 5, height = 4, unit = "in", dpi = 300)
  }
}

##Trying to make better function 3/20/21
##Time_units_factor is a factor to divide the slope by in order to make it in second units (ex if timepoints are in minutes, that time_units_factor needs to be 60)
#abort is true/false to indicate whether my_data includes abortive ligation and calculates V/E accordingly
Graphing_Kinetic_Experiments <- function(my_data, x, y, adj_y, points_to_remove, time_units_factor, graph_title, y_title, x_title, num_reps, y_lim, y_by, x_lim, x_by, base_file_name, base_folder_name, substrate, substrate_name, avg_groupings, exp, data_frame_name = "", replace_files = FALSE, filter = 0.5, z_threshold = 3, r_threshold = 0.985, V_over_E_file_name = base_file_name, y_low = 0, rm_outliers = T){
  #Globalize the arguments so they are present for functions called within this function
  for (i in arg_list){
    globalize_arguments(i, get(i))
  }
  
  ##Set Shape variable for point shapes of each rep
  shape_num <- c(15,17,18,25,3,7,8,9,10,11,12,13,14,5,2)
  my_data$shape <- NA
  for (i in 1:num_reps){
    my_data$shape[which(my_data$replicate == i)] <- shape_num[i]
  }
  
  #Removing fraction product points above filter value (default is 0.5) when the x variable is time points (not for AST)
  my_data[,adj_y] <- my_data[,y]
  my_data[points_to_remove,adj_y] <- NA
  if (x == "time_points"){
    my_data[which(my_data[,adj_y] > filter),adj_y] <- NA
  }
  my_data$removed_factor <- factor(is.na(my_data[,adj_y]), levels = c(F, T))  
  
  
  if (exp == "AST"){ #When the data is AST data
    ##Graphing individual replications
    
    #Initiate breakpoint, endpoint, rsquared values and outlier values (the return of graphing_one_rep_ast)
    my_data$breakpoint <- NA
    my_data$endpoint <- NA
    my_data$rsquared <- NA
    all_outliers <- list()
    
    #Loop through each replicate to graph it using graphing_one_rep_ast
    for(rep in unique(my_data$replicate)){
      results <- graphing_one_rep_ast(rep, my_data)
      #Grab the individual results from the return and add them to the main dataset
      all_outliers[[rep]] <- results[[1]]
      calculated_values <- results[[2]]
      my_data$breakpoint[which(my_data$replicate == rep)[1]] <- calculated_values[1]
      my_data$endpoint[which(my_data$replicate == rep)[1]] <- calculated_values[2]
      my_data$rsquared[which(my_data$replicate == rep)[1]] <- calculated_values[3]
      if (origin_manual == TRUE){
        my_data[which(my_data$replicate == rep & my_data[,x] == 0),y] <- 0
      }
    }
    #outlier_index is a list of the indices for separate datasets of each replicate
    #corr_index adds the appropriate value to convert it to the correct index for the main dataset (my_data)
    corr_index <- 0
    if (length(all_outliers) > 0) {
      for (rep in 1:length(all_outliers)){
        all_outliers[[rep]] <- all_outliers[[rep]] + corr_index
        corr_index <- corr_index + length(which(my_data$replicate == rep))
      }
    }

    ##Graphing Averages of Replicates
    
    #Initialize group variable to average across and the removed factor based on outliers from graph_one_rep_ast
    my_data$group <- paste0(my_data$time_points, "_", my_data[,x])
    my_groups <- unique(my_data$group)
    my_data$removed_factor[unlist(all_outliers)] <- as.factor(T)
    my_data[unlist(all_outliers),adj_y] <- NA
    
    ##Filtering endpoints that are lower than the filter arguement
    endpoints <- na.omit(unique(my_data$endpoint))
    low_end_reps <- c()
    if(any(endpoints < filter)){
      low_endpoints <- endpoints[which(endpoints < filter)]
      for (end in low_endpoints) {
        low_end_reps <- c(low_end_reps, my_data$replicate[which(my_data$endpoint == end)])
      }
    }
    #removing the replicates that have a low enough endpoint
    rm_reps_from_avg <- low_end_reps
    all_reps <- unique(my_data$replicate)
    new_rep_group <- all_reps
    if (length(rm_reps_from_avg) > 0){
      for (rep in rm_reps_from_avg){
        new_rep_group <- new_rep_group[new_rep_group != rep]
      }
      avg_groupings <- append(avg_groupings, new_rep_group)
    }
    #looping through each grouping of replicates being averaged
    for(grouping in 1:length(avg_groupings)){
      indices <- c()
      points_to_remove <- list_of_points_to_remove[[grouping]]
      #grabbing the indices for the replicates included
      for(w in 1:length(avg_groupings[[grouping]])){
        indices <- c(indices, which(my_data$replicate == avg_groupings[[grouping]][w]))
      }
      #Defining graphing_data based on above indices of my_data
      graphing_data <- my_data[indices,]
      
      #Initialize average and sd variables and removed for average factor and shape
      graphing_data$avg_prod_frac <- NA
      graphing_data$sd_prod_frac <- NA
      
      graphing_data$shape_for_avg <- NA
      graphing_data$removed_factor_for_avg <- factor(F, levels = c(F, T))
      
      #loop through the 'my_groups' to calculate the averages
      for (i in my_groups){
        graphing_data <- calculate_averages(i, graphing_data, exp)
      }
      
      #Create the output dataset for this set of replicates
      reps <- paste0(avg_groupings[[grouping]], collapse = "_and_")
      assign(paste0(data_frame_name, "_reps_", reps), graphing_data, envir = .GlobalEnv)
      
      #Graph this set of averages
      graphing_averages_ast(rep = avg_groupings[[grouping]], graphing_data, exp, grouping)
    } 
  }
  if (exp == "MTO_sub_dep"){
    ##Set up my_data column with the groups to calculate slopes
    my_data$slope <- NA
    my_data$V_over_E <- NA
    my_groups <<- unique(paste0("rep_", my_data$replicate, "_", substrate_name, "_", my_data[,substrate]))
    for (i in my_groups){
      my_data <- calculate_slopes(i, my_data, exp)
    }
    

    ##Graph each rep Time course
    reps_substrate_Conc <<- unique(paste0(my_data$replicate, "_", my_data[,substrate]))
    invisible(lapply(1:max(num_reps), graph_one_rep, my_data, num_reps))
    
    #Initialize appropriate variables
    my_data$adj_V_over_E <- my_data$V_over_E
    my_data$shape_V <- 16
    my_data$removed_factor_V <- factor(F, levels = c(F,T))
    my_data$kcat <- NA
    my_data$Km <- NA
    
    #Graph the Michaelis-Menten curve for each replicate
    invisible(lapply(1:num_reps, graph_one_V_over_E, my_data, exp))
    
    ##Calculate averages
    
    #Create groups based on what needs to be averaged together
    my_data$group <- paste0(my_data[,x], "_", my_data[,substrate])
    my_groups <- paste0("0_", unique(my_data[,substrate]))
    
    #See overlapping comments for if(exp == "AST") above
    if (length(avg_groupings) > 0){
      for (grouping in 1:length(avg_groupings)){
        indices <- c()
        points_to_remove <- list_of_points_to_remove[[grouping]]
        for(w in 1:length(avg_groupings[[grouping]])){
          indices <- c(indices, which(my_data$replicate == avg_groupings[[grouping]][w]))
        }
        graphing_data <- my_data[indices,]

        graphing_data$avg_V_over_E <- NA
        graphing_data$sd_V_over_E <- NA
        
        graphing_data$adj_V_over_E_for_avg <- graphing_data$V_over_E
        graphing_data$shape_V_for_avg <- 16
        graphing_data$removed_factor_V_for_avg <- factor(F, levels = c(F,T))
        
        for (i in my_groups){
          graphing_data <- calculate_averages(i, graphing_data, exp)
        }
        
        ##Graph the V/E averages
        graphing_data$removed_factor_avg <- factor(F, levels = c(F, T))
        graphing_data$avg_adj_V_over_E <- graphing_data$avg_V_over_E
        graphing_data$shape_avg <- 16
        
        graph_averages_V_over_E(avg_groupings[[grouping]], graphing_data, exp)
        
        ##Save the dataset for this set of replicates
        print(avg_groupings[[grouping]])
        reps <- paste0(avg_groupings[[grouping]], collapse = "_and_")
        assign(paste0(data_frame_name, "_reps_", reps), graphing_data, envir = .GlobalEnv)
      }
    }
  }

  #remove the globalized arguments so they don't mess up future runs
  for (i in arg_list){
    remove_globalized_arguments(get(i))
  }
}


#Organizing Data Function

reading_in_data <- function(filename, date_index = 1)
  #filename: Name of the file being read in
  #date_index: If filenames include the date, separate it from the rest of the name with spaces. The date_index indicates the folder level the files are in from the working directory
  {
  #Read in the file: Assumes tab delimited and with a header (format from the ImageQuant text files)
  data <- as.data.frame(t(read.table(filename, header = T, stringsAsFactors = F, sep = "\t")))
  #First column (after transpose) is the 'Lane number' so remove
  data$V1 <- NULL
  #For example data the first band is product and the second band is substrate
  colnames(data) <- c("product", "substrate")
  
  #Ensure everything is being treated as a numeric
  data <- as.data.frame(apply(data, 2, as.numeric))
  
  #calculate the fraction product value
  data$prod_frac <- data$product/(data$product + data$substrate)
  
  #Grab the date from the file name and add a columnin the dataset
  date <- strsplit(strsplit(filename, split = "/")[[1]][date_index], split = " ")[[1]][1]
  data$date_done <- date
  
  #Ongoing table of the dimensions for each file to use when determining indices_wanted
  dim_table <<- rbind(dim_table, data.frame(file_rows = dim(data)[1], cuml_rows = (dim(data)[1]) + sum(dim_table[,1]), files_name = filename))
  #returns the data
  return(data)
}

#Essentially takes uses the organization of the arguments to populate the dataset with the experimental conditions
#See Template Script.Rmd for details on how to organize the arguments for an example dataset
organizing_data <- function(data, indices_wanted, experiment, proteins, num_reps, num_rxns_per_rep, num_time_pts_each, time_points, Enz_Conc, DNA_Conc, Mg_Conc, ATP_Conc, substrates = c(0), mult = FALSE)
  #data is the output from organize_data() 
  #indices_wanted is the used to order the data into the order the rest of the variables are in
  #experiment is the name of the experiment
  #proteins a vector of the proteins tested
  #num_reps is a vector or list for the number of reps for each experimental set
  #num_rxns_per_rep is a vector or list for the number of reactions per replicate
  #num_time_pts_each is a vector or list for the number of time points per reaction
  #time_points is a list of the time points taken for each replicate
  #Enz_Conc is a list of the enzyme concentration for each replicate
  #DNA_Conc is a list of the DNA concentration for each replicate 
  #Mg_Conc is a list of the Mg concentration for each replicate
  #ATP_Conc is a list of the ATP concentration for each replicate
  #substrates is a vector of the substrates used (default is c())
  #mult is whether there are multiple substrates (default is FALSE)
  {
  data <- data[indices_wanted,]
  
  num_prot <- length(proteins)
  num_sub <- length(substrates)
  if (mult){
    num_protein_each <- c()
    substrate_var <- c()
    replicate_vector <- c()
    times <- c()
    DNAs <- c()
    Enzs <- c()
    ATPs <- c()
    Mgs <- c()
    for (rep in 1:max(unlist(num_reps))){
      for (sub in 1:num_sub){
        for (prot in 1:num_prot){
          num_protein_each <- c(num_protein_each, num_rxns_per_rep[[sub]][[rep]][prot]*num_time_pts_each[[sub]][[rep]][prot])
          substrate_var <- c(substrate_var, rep(substrates[sub], num_rxns_per_rep[[sub]][[rep]][prot]*num_time_pts_each[[sub]][[rep]][prot]))
          replicate_vector <- c(replicate_vector, rep(rep, num_rxns_per_rep[[sub]][[rep]][prot]*num_time_pts_each[[sub]][[rep]][prot]))
          times <- c(times, time_points[[rep]][[sub]][[prot]])
          
          each <- num_time_pts_each[[sub]][[rep]][prot]
          DNAs <- c(DNAs, rep(DNA_Conc[[rep]][[sub]][[prot]], each = each))
          Enzs <- c(Enzs, rep(Enz_Conc[[rep]][[sub]][[prot]], each = each))
          ATPs <- c(ATPs, rep(ATP_Conc[[rep]][[sub]][[prot]], each = each))
          Mgs <- c(Mgs, rep(Mg_Conc[[rep]][[sub]][[prot]], each = each))
        }
      }
    }
    
    data$protein <- NA
    for (i in 1:length(num_protein_each)){
      if (num_protein_each[i]==0){
        next
      } else {
        start_index <- sum(num_protein_each[0:(i-1)])+1
        end_index <- sum(num_protein_each[1:i])
        num_protein <- length(proteins)
        if (i%%num_protein == 0){
          data$protein[start_index:end_index] <- proteins[num_protein]
        } else{
          data$protein[start_index:end_index] <- proteins[i%%num_protein]
        }
      }
    }
    
    data$Substrate <- substrate_var
    data$replicate <- replicate_vector
    data$time_points <- times
    data$DNA_Conc <- DNAs
    data$Enz_Conc <- Enzs
    data$ATP_Conc <- ATPs
    data$Mg_Conc <- Mgs
    
    if (grepl("MTO", experiment) | experiment == "AST_no_control"){
      #browser()
      beg_point_each <- c() 
      for (rep in 1:max(unlist(num_reps))){
        for (sub in 1:num_sub){
          for (prot in 1:num_prot){
            if (rep == 1 & sub == 1 & prot == 1){
              start_point <- 1
              end_point <- num_rxns_per_rep[[sub]][[rep]][prot] * (num_time_pts_each[[sub]][[rep]][prot]+1)
              beg_point_each <- c(beg_point_each, seq(from = start_point, to = end_point, by = (num_time_pts_each[[sub]][[rep]][prot]+1)))
              new_start_point <- end_point + 1
            } else {
              start_point <- new_start_point
              end_point_new <- num_rxns_per_rep[[sub]][[rep]][prot] * (num_time_pts_each[[sub]][[rep]][prot]+1) + end_point
              if (start_point >= end_point_new){
                next
              } else {
                end_point <- end_point_new
                beg_point_each <- c(beg_point_each, seq(from = start_point, to = end_point, by = (num_time_pts_each[[sub]][[rep]][prot]+1)))
                new_start_point <- end_point + 1
              }
            }
          }
        }
      }
      
      #browser()  
      for (i in beg_point_each){
        num_col <- ncol(data)
        num_zeros <- which(colnames(data) == "date_done") - 1
        data <- InsertRow(data, c(rep(0,num_zeros),data[i,(num_zeros+1):num_col]), i)
        data$time_points[i] <- 0
      }
    }
    
  } else {
    num_protein_each <- c()
    replicate_vector <- c()
    times <- c()
    DNAs <- c()
    Enzs <- c()
    ATPs <- c()
    Mgs <- c()
    for (rep in 1:max(unlist(num_reps))){
      for (prot in 1:num_prot){
        num_protein_each <- c(num_protein_each, num_rxns_per_rep[[rep]][prot]*num_time_pts_each[[rep]][prot])
        replicate_vector <- c(replicate_vector, rep(rep, num_rxns_per_rep[[rep]][prot]*num_time_pts_each[[rep]][prot]))
        times <- c(times, time_points[[rep]][[prot]])
          
        each <- num_time_pts_each[[rep]][prot]
        DNAs <- c(DNAs, rep(DNA_Conc[[rep]][[prot]], each = each))
        Enzs <- c(Enzs, rep(Enz_Conc[[rep]][[prot]], each = each))
        ATPs <- c(ATPs, rep(ATP_Conc[[rep]][[prot]], each = each))
        Mgs <- c(Mgs, rep(Mg_Conc[[rep]][[prot]], each = each))
      }
    }
    
    data$protein <- NA
    for (i in 1:length(num_protein_each)){
      if (num_protein_each[i]==0){
        next
      } else {
        start_index <- sum(num_protein_each[0:(i-1)])+1
        end_index <- sum(num_protein_each[1:i])
        num_protein <- length(proteins)
        if (i%%num_protein == 0){
          data$protein[start_index:end_index] <- proteins[num_protein]
        } else{
          data$protein[start_index:end_index] <- proteins[i%%num_protein]
        }
      }
    }
    
    data$replicate <- replicate_vector
    data$time_points <- times
    data$DNA_Conc <- DNAs
    data$Enz_Conc <- Enzs
    data$ATP_Conc <- ATPs
    data$Mg_Conc <- Mgs
    
    if (grepl("MTO", experiment) | experiment == "AST_no_control"){
      #browser()
      beg_point_each <- c() 
      for (rep in 1:max(unlist(num_reps))){
        for (prot in 1:num_prot){
          if (rep == 1 & sub == 1 & prot == 1){
            start_point <- 1
            end_point <- num_rxns_per_rep[[rep]][prot] * (num_time_pts_each[[rep]][prot]+1)
            beg_point_each <- c(beg_point_each, seq(from = start_point, to = end_point, by = (num_time_pts_each[[rep]][prot]+1)))
            new_start_point <- end_point + 1
          } else {
            start_point <- new_start_point
            end_point_new <- num_rxns_per_rep[[rep]][prot] * (num_time_pts_each[[rep]][prot]+1) + end_point
            if (start_point >= end_point_new){
              next
            } else {
              end_point <- end_point_new
              beg_point_each <- c(beg_point_each, seq(from = start_point, to = end_point, by = (num_time_pts_each[[rep]][prot]+1)))
              new_start_point <- end_point + 1
            }
          }
        }
      }
      
      #browser()  
      for (i in beg_point_each){
        num_col <- ncol(data)
        num_zeros <- which(colnames(data) == "date_done") - 1
        data <- InsertRow(data, c(rep(0,num_zeros),data[i,(num_zeros+1):num_col]), i)
        data$time_points[i] <- 0
      }
    }
  }

  data$prod_frac <- as.numeric(data$prod_frac)
  data$replicate <- as.numeric(data$replicate)
  data$time_points <- as.numeric(data$time_points)
  return(data)
}




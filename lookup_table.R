#-----------------------
#     LOAD PACKAGES
#-----------------------
library(ncdf4)
library(reshape2)

#-----------------------
#     Route
#-----------------------
met_route <- '/home/eric/era-interim/for_mtheta/'
table_nh_route <- '/home/eric/yuming/m_theta/era_nh_pv_2/'
table_sh_route <- '/home/eric/yuming/m_theta/era_sh_pv_2/'
table_route <- '/home/eric/yuming/m_theta/era_lookup_merge/'

#-----------------------
#    Define Global Parameter
#-----------------------

authalic_radius_m <- (6371.0072 * 1000) # radius of earth in m
Mw_air_g_mol <- 28.97 # molecular weight of dry air in g/mol      
R_mbar_m3_K_mol <- (8.3144598e-5)*1000 # ideal gas constant in mbar m3 / K mol

Rd <- 287.04 #specific gas constant for air (287.04 J/(kg·K))
Cpd <- 1005.7 #specific heat of dry air at constant pressure (1005.7 J/(kg·K)).

pv_cut <- 2e-6 # Define tropopause


#-----------------------
#     Compute Mthetae
#-----------------------
for (yr in 1980:2018) {
  setwd(met_route)
  nc <- nc_open(paste(yr,'.nc',sep=''))
  date <- seq(as.Date(paste(yr,'01','01',sep='-')),as.Date(paste(yr,'12','31',sep='-')),by="day")
  theta_input <- seq(220,420,1)
  lookup <- data.frame(matrix(NA,nrow=length(date),ncol = length(theta_input)+1))
  colnames(lookup) <- c('date',theta_input)
  lookup$date <- date
  lookup_nh <- lookup
  lookup_sh <- lookup
  rm(lookup)
  lon <- ncvar_get(nc=nc,varid='longitude')
  lon <- lon - 180 #convert to -180 to 177.5
  lat <- ncvar_get(nc=nc,varid='latitude')
  lat <- rev(lat) # convert to -90 to 90
  p <- ncvar_get(nc=nc,varid = 'level')
  p <- rev(p)
  sfp <- nc_open(paste('sfp',yr,'.nc',sep=''))
  for (dy in 1:length(date)) {
    
    #Read Portential Vorticity, Relative Humidity, Temperature and surface pressure mask
    pv <- ncvar_get(nc=nc,varid='pv',start = c(1,1,1,(4*dy-3)),count = c(length(lon),length(lat),length(p),4))
    pv <- apply(pv,c(1,2,3),mean)
    pv <- pv[c(73:144,1:72),73:1,30:1]
    temp <- ncvar_get(nc=nc,varid='t',start = c(1,1,1,(4*dy-3)),count = c(length(lon),length(lat),length(p),4))    
    temp <- apply(temp,c(1,2,3),mean)
    temp <- temp[c(73:144,1:72),73:1,30:1]
    rh <- ncvar_get(nc=nc,varid='r',start = c(1,1,1,(4*dy-3)),count = c(length(lon),length(lat),length(p),4))    
    rh <- apply(rh,c(1,2,3),mean)
    rh <- rh[c(73:144,1:72),73:1,30:1]
    sp <- ncvar_get(nc=sfp,varid='sp',start = c(1,1,(4*dy-3)),count = c(length(lon),length(lat),4))    
    sp <- apply(sp,c(1,2),mean)
    sp <- sp[c(73:144,1:72),73:1]/100
    
    #Compute mass
    half_lon <- (unique(diff(lon))/2)[1]
    half_lat <- (unique(diff(lat))/2)[1]
    
    ### Compute area of each grid cell (2D)
    lat0 <- lat-half_lat
    #lat0[1] <- -90
    lat0 <- lat0*(pi/180)
    lat1 <- lat+half_lat
    #lat1[length(lat1)] <- 90
    lat1 <- lat1*(pi/180)
    
    lon0 <- (lon-half_lon)*(pi/180)
    lon1 <- (lon+half_lon)*(pi/180)
    
    area_m2 <- array(data=NA, dim=c(length(lon), length(lat))) # create empty array
    
    for(i in seq(from=1, to=length(lon))){
      for(j in seq(from=1, to=length(lat))){
        
        area_m2[i,j] <- abs(sin(lat1[j]) - sin(lat0[j])) *
          abs(lon1[i] - lon0[i]) * (authalic_radius_m^2)
        
      }
    }
    
    ### Calculate mass per m2 at each pressure height section (do not need to have height)
    p <- rev(p)
    half_p <- c(-p[1],p[1:(length(p)-1)] - p[2:length(p)])/2
    tops_p <- p+half_p
    cube_p_diff <- c(diff(tops_p), -tops_p[length(tops_p)]+1013.25)
    p <- rev(p)
    cube_p_diff <- rev(cube_p_diff)
    tops_p <- rev(tops_p)
    
    ###Create area , pressure and pressure difference cube
    area_m2_3d <- array(data=NA, dim=c(length(lon), length(lat), length(p)))
    
    for(q in seq(from=1, to=length(cube_p_diff))){
      
      area_m2_3d[,,q] <- area_m2
      
    }
    
    cube_p_diff_3d <- array(data=NA, dim=c(length(lon), length(lat), length(p)))
    
    for(q in seq(from=1, to=length(cube_p_diff))){
      
      cube_p_diff_3d[,,q] <- cube_p_diff[q]
      
    }
    
    p_mbar_3d <- array(data=NA, dim=c(length(lon), length(lat), length(p)))
    
    for(q in seq(from=1, to=length(cube_p_diff))){
      
      p_mbar_3d[,,q] <- p[q]
      
    }
    
    # Create g cube 
    g_array <- array(data=NA,dim=c(length(lon), length(lat), length(p)))
    for( q in 1:length(lat)){
      g_array[,q,] <- lat[q] * pi / 180
    }
    g_array <- 9.78046 * (1+0.0052884*(sin(g_array))^2-0.0000059*(sin(2*g_array))^2)
    # height_m <- ((1 - ((p / 1013.25) ^ 0.190284)) * 145366.45) * 0.3048
    height_m <- log(p/1013.25)*(-8400)
    for ( q in 1:length(height_m)) {
      g_array[,,q] <- g_array[,,q] - 0.000003086*height_m[q]
      
    }
    
    ### now compute mass of each grid cell in g. mass
    mass_g <- cube_p_diff_3d * 100 *area_m2_3d / g_array
    
    
    ### Calculate saturation mixing ratio
    #Begin with saturation vapor pressure
    p_sv <- 6.112*exp(17.67*(temp-273.15)/((temp-273.15)+243.5))
    
    #saturation mixing ratio
    ws <- 0.622*p_sv/(p_mbar_3d)
    
    #mixing ratio
    w <- rh*ws/100
    
    #mass of water vapor
    m_water <- mass_g * w
    
    #mass of dry air
    mass <- mass_g - m_water
    
    #compute thetae  
    theta_e <- array(data=NA, dim=dim(temp)) 
    for(i in 1:length(p)){
      theta_e[,,i] <- (temp[,,i]+ (array(data=2501,dim=dim(sp))-((2501-2406)/40)*(temp[,,i]-273.15))*1000*(w[,,i]/Cpd))*(1013/p[i])^(Rd/Cpd)
    }
    
    #---------------------------------------
    #     Assgian Dimnames
    #---------------------------------------
    dimnames(theta_e) <- list(lon,lat,p)
    dimnames(mass) <- list(lon,lat,p)
    dimnames(pv) <- list(lon,lat,p)
    
    #---------------------------------------
    #     Cut off by using PVU=2
    #---------------------------------------
    
    for (lon_index in 1:length(lon)) {
      for (lat_index in 1:length(lat)) {
        p_range <- which(p <= 400)
        
        p_index <- min(which(pv[lon_index,lat_index,p_range] > pv_cut | pv[lon_index,lat_index,p_range] < -pv_cut))
        if(!is.infinite(p_index)){
          mass[lon_index,lat_index,p_range[p_index]:length(p)] <- NA 
        }
      }
    }
    
    #---------------------------------------
    #     Consider surface topography
    #---------------------------------------
    for (i in 1:length(lon)) {
      for (j in 1:length(lat)) {
        change_index <- which(tops_p > sp[i,j])
        if(length(change_index) >=  1){
          theta_e[i,j,change_index] <- NA
          pv[i,j,change_index] <- NA
          mass[i,j,change_index] <- NA
          mass[i,j,max(change_index)+1] <- mass[i,j,max(change_index)+1]* (sp[i,j] - tops_p[max(change_index)+1]) /
            (tops_p[max(change_index)] - tops_p[max(change_index)+1] )
        }else{
          mass[i,j,1] <- mass[i,j,1] * (sp[i,j]-tops_p[1]) / (1013.25 - tops_p[1])
        }
      }
    }
    
    check_table <- array(1,dim=dim(pv))
    check_table <- check_table*mass
    check_table[!is.na(check_table)] <- 1
    theta_e <- theta_e*check_table
    
    #---------------------------------------
    #     Divide into NH & SH
    #---------------------------------------
    
    theta_e_north <- theta_e[,which(lat > 0),]
    
    mass_north <- mass[,which(lat > 0),]
    
    pv_north <- pv[,which(lat > 0),]
    
    theta_e_south <- theta_e[,which(lat < 0),]
    
    mass_south <- mass[,which(lat < 0),]
    
    pv_south <- pv[,which(lat < 0),]
    
    
    #---------------------------------------
    #     SORT BY THETA/CUMULATIVE MASS Create lookuptable
    #---------------------------------------
   #Northern Hemisphere
    setwd(table_nh_route) 
    
    df <- melt(mass_north)
    tmp <- melt(theta_e_north)
    df <- cbind(df,tmp[,4])
    tmp <- melt(pv_north)
    df <- cbind(df,tmp[,4])
    colnames(df) <- c("longitude","latitude","pressure","mass","theta_e",'pv')
    df <- df[order(df$theta_e),]
    df$m_theta_e <- cumsum(df$mass)
    
    write.csv(df,file = paste(date[dy],'_north.csv',sep=''))
    m_theta_e <- approx(x=df$theta_e,y=df$m_theta_e,xout=theta_input,rule=1,method = 'linear')$y
    lookup_nh[dy,2:ncol(lookup_nh)] <- m_theta_e
    
    
    #Southern Hemisphere
    setwd(table_sh_route) 
    
    df <- melt(mass_south)
    tmp <- melt(theta_e_south)
    df <- cbind(df,tmp[,4])
    tmp <- melt(pv_south)
    df <- cbind(df,tmp[,4])
    colnames(df) <- c("longitude","latitude","pressure","mass","theta_e",'pv')
    df <- df[order(df$theta_e),]
    df$m_theta_e <- cumsum(df$mass)
    
    write.csv(df,file = paste(date[dy],'_south.csv',sep=''))
    m_theta_e <- approx(x=df$theta_e,y=df$m_theta_e,xout=theta_input,rule=1,method = 'linear')$y
    lookup_sh[dy,2:ncol(lookup_sh)] <- m_theta_e
    print(paste(date[dy],'-finised',sep=''))
  }
  
  #save annual merge table
  setwd(table_route)
  write.csv(lookup_nh,file = paste(yr,'_north.csv',sep=''))
  write.csv(lookup_sh,file = paste(yr,'_south.csv',sep=''))
}

#merge table
for (year in 1980:2018) {
  lookup <- read.csv(paste('/home/eric/yuming/m_theta/era_lookup_merge/',year,'_north.csv',sep=''))
  lookup <- lookup[,-1]
  if(year == 1980){
    lookup_merge <- lookup
  }else{
    lookup_merge <- rbind(lookup_merge,lookup)
  }
}

write.table(lookup_merge,'/home/eric/yuming/m_theta/era_lookup_merge/lookuptable_NH.txt',row.names = F)


for (year in 1980:2018) {
  lookup <- read.csv(paste('/home/eric/yuming/m_theta/era_lookup_merge/',year,'_south.csv',sep=''))
  lookup <- lookup[,-1]
  if(year == 1980){
    lookup_merge <- lookup
  }else{
    lookup_merge <- rbind(lookup_merge,lookup)
  }
}

write.table(lookup_merge,'/home/eric/yuming/m_theta/era_lookup_merge/lookuptable_SH.txt',row.names = F)
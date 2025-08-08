```R
mod.nls.nostart <- function(data, TFoi){
  mod.nls <- nls(Overall_transcriptomic_change ~ SSlogis(Dose,  Asym, xmid, scal), 
                 data = data, trace = F)
  return(mod.nls)
}

mod.nls.start.range <- function(data, TFoi, init, upper.value, lower.value){
  mod.nls <- nls(Overall_transcriptomic_change ~ SSlogis(Dose,  Asym, xmid, scal), 
                 data = data,
                #  data = subset(df, TF == TFoi),
                 start = init,
                 upper = upper.value,
                 lower = lower.value,
                 algorithm = "port", trace = F)
  return(mod.nls)
}

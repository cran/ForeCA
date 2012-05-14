discrete.entropy <-
function(prob, base = exp(1), method = "MLE", threshold = 0){
  
  prob = prob / sum(prob)
  
  if (threshold > 0) {
    prob[prob < threshold] = 0
  }
  if (method == "smoothed"){
    prob = prob + 1/length(prob)
    prob = prob / sum(prob)
    HH = - sum(prob * log(prob, base = base))
  }
  if (method == "MLE"){
    if (any(prob == 0)) {
      prob = prob[prob != 0]
    }
    HH = -sum(prob * log(prob, base = base))
  }

  attr(HH, "base") = as.character(base)
  return(HH)
}

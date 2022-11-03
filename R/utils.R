#'@title MSE
#'@export
mse = function(x,y){mean((x-y)^2)}

#'@title RMSE
#'@export
rmse = function(x,y){sqrt(mean((x-y)^2))}

#'@title MAE
#'@export
mae = function(x,y){mean(abs(x-y))}

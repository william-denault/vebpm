#'@title MSE
mse = function(x,y){mean((x-y)^2)}

#'@title RMSE
rmse = function(x,y){sqrt(mean((x-y)^2))}

#'@title MAE
mae = function(x,y){mean(abs(x-y))}

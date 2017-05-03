function dudt = lecture7_noncon_burgers_ODE(t,u,D)
  dudt = -u.*(D*u);
end

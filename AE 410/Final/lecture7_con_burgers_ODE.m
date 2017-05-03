function dudt = lecture7_con_burgers_ODE(t,u,D)
  dudt = -0.5*D*(u.*u);
end

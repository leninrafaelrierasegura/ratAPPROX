

matern.p.precision(loc = c(0,1), 
                   kappa = 10, 
                   p = 0,
                   equally_spaced = TRUE, 
                   alpha = 1.8)$Q

r00_inverse <- matern.p.joint(s = 0, t = 0, kappa = 10, p = 0, alpha = 1.8)
r01_inverse <- matern.p.joint(s = 0, t = 1, kappa = 10, p = 0, alpha = 1.8)
r10_inverse <- matern.p.joint(s = 1, t = 0, kappa = 10, p = 0, alpha = 1.8)
r11_inverse <- matern.p.joint(s = 1, t = 1, kappa = 10, p = 0, alpha = 1.8)

solve(rbind(cbind(r00_inverse, r01_inverse),
            cbind(r10_inverse, r11_inverse)))
xi = 0.57735
eta = -0.57735
X = [8, -0.5; 8.5, -0.5; 8.5, 0.; 8, 0; 8.25, -0.5; 8.5, -0.25; 8.25, 0.; 8, -0.25]
true = 1

[B, detJ] = B_Quad8(xi, eta, X, true)

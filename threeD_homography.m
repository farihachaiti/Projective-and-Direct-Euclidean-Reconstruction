function H=threeD_homography(X1p2,X2p2,X3p2,X4p2,X5p2,X1E,X2E,X3E,X4E,X5E)
  clc

   #Measuring corresponding points
  printf("Measuring corresponding points\n")
 
  [X1p2c, X2p2c, X3p2c, X4p2c,X5p2c,T1] = measure_points(X1p2,X2p2,X3p2,X4p2,X5p2,0)
[X1Ec,X2Ec,X3Ec,X4Ec,X5Ec,T2] = measure_points(X1E,X2E,X3E,X4E,X5E,1)

  #Calculating Design Matrix
  printf("Calculating Design Matrix\n")
 
 A1 = design_matrix(X1p2c,X1Ec)
 A2 = design_matrix(X2p2c,X2Ec)
 A3 = design_matrix(X3p2c,X3Ec)
 A4 = design_matrix(X4p2c,X4Ec)
 A5 = design_matrix(X5p2c,X5Ec)
 
 A = [A1; A2; A3 ;A4 ;A5 ;0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
 
  #Calculating 3D Homography
  printf("Calculating 3D Homography\n")
 
 
 [U,D,V] = svd(A)
 h = [V(:,end)]
 
  Hconditioned = transpose(reshape(h,[4,4]))
  
  H = inverse(T2)*Hconditioned*T1
  
 
 end


#function [y,c]=define_function(x)
function [Y1,Y2,Y3,Y4,Y5,T]=measure_points(X1,X2,X3,X4,X5,E)
  
   #Conditioning: Translation
  printf("Conditioning: Translation\n")
  
  XYZmat = [abs(X1) abs(X2) abs(X3) abs(X4) abs(X5)]
  
  t = mean(XYZmat,2)
  t(1,1)
  t(2,1)
  t(3,1)
  
   #Conditioning: Scalling
  printf("Conditioning: Scalling\n")
  
 X1xfinal = X1(1,1) - t(1,1)
 X1yfinal = X1(2,1) - t(2,1)
 X1zfinal = X1(3,1) - t(3,1)
 X2xfinal = X2(1,1) - t(1,1)
 X2yfinal = X2(2,1) - t(2,1)
 X2zfinal = X2(3,1) - t(3,1)
 X3xfinal = X3(1,1) - t(1,1)
 X3yfinal = X3(2,1) - t(2,1)
 X3zfinal = X3(3,1) - t(3,1)
 X4xfinal = X4(1,1) - t(1,1)
 X4yfinal = X4(2,1) - t(2,1)
 X4zfinal = X4(3,1) - t(3,1)
 X5xfinal = X5(1,1) - t(1,1)
 X5yfinal = X5(2,1) - t(2,1)
 X5zfinal = X5(3,1) - t(3,1)
  
 if(E==0)
 X1Y1Z1final = [X1xfinal; X1yfinal; X1zfinal; X1(4,1) ]
 X2Y2Z2final = [X2xfinal; X2yfinal; X2zfinal; X2(4,1) ]
 X3Y3Z3final = [X3xfinal; X3yfinal; X3zfinal; X3(4,1) ]
 X4Y4Z4final = [X4xfinal; X4yfinal; X4zfinal; X4(4,1) ]
 X5Y5Z5final = [X5xfinal; X5yfinal; X5zfinal; X5(4,1) ]
else
 X1Y1Z1final = [X1xfinal; X1yfinal; X1zfinal; 1 ]
 X2Y2Z2final = [X2xfinal; X2yfinal; X2zfinal; 1 ]
 X3Y3Z3final = [X3xfinal; X3yfinal; X3zfinal; 1 ]
 X4Y4Z4final = [X4xfinal; X4yfinal; X4zfinal; 1 ]
 X5Y5Z5final = [X5xfinal; X5yfinal; X5zfinal; 1 ]
endif

 XYZmat2 = [abs(X1Y1Z1final) abs(X2Y2Z2final) abs(X3Y3Z3final) abs(X4Y4Z4final) abs(X5Y5Z5final)]
 
 s = mean(XYZmat2,2)
 s(1,1)
 s(2,1)
 s(3,1)
 
  #Coordinate Transformation
  printf("Coordinate Transformation\n")
 
 
 T = [1/s(1,1) 0 0 0; 0 1/s(2,1) 0 0; 0 0 1/s(3,1) 0;0 0 0 1]*[1 0 0 -t(1,1); 0 1 0 -t(2,1); 0 0 1 -t(3,1); 0 0 0 1] 
 
  #Conditioned Coordinates
  printf("Conditioned Coordinates\n")
 
 
 Y1 = T*X1Y1Z1final
 Y2 = T*X2Y2Z2final
 Y3 = T*X3Y3Z3final
 Y4 = T*X4Y4Z4final
 Y5 = T*X5Y5Z5final
 
endfunction

function y=design_matrix(Xp2c,XEc)
   y = [-XEc(4,1)*Xp2c(1,1)  -XEc(4,1)*Xp2c(2,1) -XEc(4,1)*Xp2c(3,1) -XEc(4,1)*Xp2c(4,1) 0 0 0 0 0 0 0 0 XEc(1,1)*Xp2c(1,1) XEc(1,1)*Xp2c(2,1) XEc(1,1)*Xp2c(3,1) XEc(1,1)*Xp2c(4,1);0 0 0 0 -XEc(4,1)*Xp2c(1,1) -XEc(4,1)*Xp2c(2,1) -XEc(4,1)*Xp2c(3,1) -XEc(4,1)*Xp2c(4,1) 0 0 0 0 XEc(2,1)*Xp2c(1,1) XEc(2,1)*Xp2c(2,1) XEc(2,1)*Xp2c(3,1) XEc(2,1)*Xp2c(4,1);0 0 0 0 0 0 0 0 -XEc(4,1)*Xp2c(1,1) -XEc(4,1)*Xp2c(2,1) -XEc(4,1)*Xp2c(3,1) -XEc(4,1)*Xp2c(4,1) XEc(3,1)*Xp2c(1,1) XEc(3,1)*Xp2c(2,1) XEc(3,1)*Xp2c(3,1) XEc(3,1)*Xp2c(4,1)]
  
  endfunction


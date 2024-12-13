function mainfunction
  clc
  
  #Taking  image points
  printf("Taking  image points\n")

  [xy,xyprime]=take_image_points_from_file()

  
    #Measuring image points
  printf("Measuring image points\n")
  

  [xy_final,T1] = measure_points(xy)
  [xy_prime_final,T2] = measure_points(xyprime)
  
    #Generating Design Matrix
  printf("Generating Design Matrix\n")
 
   for k = 1:length(xy)
    A{k} = design_matrix(xy_final{k},xy_prime_final{k})
   end
 
   #Calculating Fundamental Matrix F
  printf("Calculating Fundamental Matrix F\n")
 
 F1 = calculate_fundamental_matrix(A,T1,T2);
 
  #Forcing Singularity Constraint on F
  printf("Forcing Singularity Constraint on F\n")
 
 F = force_rank2(F1)
 
     #Calculating epipoles from F
  printf("Calculating epipoles from Matrix F\n")
  
  [e,ePrime]=calculate_epipoles(F)
  
  
  #Calculating Projection Matrices Pn and Pprime from F
  printf("Calculating Projection Matrices Pn and Pprime from F\n")
  
  
[Pn,Pprime]=get_projection_matrix(F,ePrime)
 
  #Calculating Linear Triangulation 
  printf("Calculating Linear Triangulation \n")
  
  for k = 1:length(xy)
    Xp1{k} = linear_triangulation(Pn,Pprime,xy{k},xyprime{k})
   end


  #Visualizing Spatial Object Points 
  printf("Visualizing Spatial Object Points \n")

  visualize_spatial_object_point(Xp1)

 

 #Taking Control Points 
  printf("Taking Control Points \n")
 
 [xy_cp,xy_cp_prime,X1E,X2E,X3E,X4E,X5E]=take_control_points_from_file()
 
  #Triangulating projective object points 
  printf("Triangulating projective object points\n")
  
   for k = 1:length(xy_cp)
    Xp2{k} = linear_triangulation(Pn,Pprime,xy_cp{k},xy_cp_prime{k})
   end


 
  #Calculating 3D Homography
  printf("Calculating 3D Homography \n")
 
 H = threeD_homography(Xp2{1},Xp2{2},Xp2{3},Xp2{4},Xp2{5},X1E,X2E,X3E,X4E,X5E)
 H = H/H(4,4)
 
 #Euclidean Reconstruction
  printf("Euclidean Reconstruction \n")
 
  for k = 1:length(Xp1)
    XEfinal{k} = H*Xp1{k}
   end
 
 for k = 1:length(Xp1)
    XEfinal2{k} = XEfinal{k}/XEfinal{k}(4,1)
   end
 
 
  #Visualizing the result of Euclidean Reconstruction spatially 
  printf("Visualizing the result of Euclidean Reconstruction spatially \n")
 

  visualize_spatial_object_point(XEfinal2)

 
 endfunction

 
 
 function visualize_spatial_object_point(X)
   figure; hold on;
   for k = 1:length(X)
    scatter3(X{k}(1,:), X{k}(2,:), X{k}(3,:), 10, 'filled')
   end
   axis square; view(32, 75);
   hold off;
   grid on;
 endfunction

 
 

 function y=design_matrix(xy,xy_prime)
   y = [xy(1,1)*xy_prime(1,1)  xy(2,1)*xy_prime(1,1) xy_prime(1,1) xy(1,1)*xy_prime(2,1) xy(2,1)*xy_prime(2,1) xy_prime(2,1) xy(1,1) xy(2,1) 1]
  
endfunction


  
 function [XY_conditioned,T]=measure_points(xy)
  
   #Conditioning: Translation
  printf("Conditioning: Translation\n")
  
 m = 3
 n = length(xy)
 XYmat = zeros(m,n);
 
 for k = 1:length(xy)
  XY{k} = [xy{k}; 1]
 end
 
 
 for i = 1:m
   for j = 1:n
    XYmat(i,j) = [abs(XY{j}(i))]
   end
 end
   
  t = mean(XYmat,2)
  t(1,1)
  t(2,1)
  t(3,1)

  
   #Conditioning: Scalling
  printf("Conditioning: Scalling\n")
  
   for k = 1:length(xy)
    Xfinal{k} = xy{k}(1,1) - t(1,1)
    Yfinal{k} = xy{k}(2,1) - t(2,1)
    XYfinal{k} = [Xfinal{k}; Yfinal{k}; 1]
   end

   XYmat2 = zeros(m,n)
 
  for i = 1:m
   for j = 1:n
    XYmat2(i,j) = [abs(XY{j}(i))]
   end
 end
   
XYmat
XYmat2
 
 s = mean(XYmat2,2)
 s(1,1)
 s(2,1)
 s(3,1)
 
  #Coordinate Transformation
  printf("Coordinate Transformation\n")
 
 
 T = [1/s(1,1) 0 0; 0 1/s(2,1) 0; 0 0 1]*[ 1 0 -t(1,1); 0 1 -t(2,1); 0 0 1] 
 
  #Conditioned Coordinates
  printf("Conditioned Coordinates\n")
 
  for k = 1:length(xy)
    XY_conditioned{k} = T* XY{k}
   end

 
endfunction




function X=linear_triangulation(P,Pprime,xy,xyprime)
   A = [xy(1,1)*P(3,:)-P(1,:); xy(2,1)*P(3,:)-P(2,:);  xyprime(1,1)*Pprime(3,:)-Pprime(1,:); xyprime(2,1)*Pprime(3,:)-Pprime(2,:)]
 
 [U,D,V] = svd(A)
 
 X = [V(:,end)]
 
endfunction




 function F=calculate_fundamental_matrix(A,T1,T2)
    
 
  #Design Matrix A1
  printf("Design Matrix A\n")
 m = length(A)
 n = 9
 Afinal = zeros(m,n)
 
 for i = 1:m
   for j = 1:n
    Afinal(i,j) = [A{i}(j)]
   end
 end
   
 
  #Singular Value Decompisition
  printf("Singular Value Decompisition\n")
 
 [U,D,V] = svd(Afinal)

 F_prime = [V(:,end)]
 
  #Reshaping to get the Fundamental Matrix
  printf("Reshaping to get the Fundamental Matrix\n")
 
 F2 = reshape(F_prime,3,3)'

   #Reverse Conditioning
  printf("Reverse Conditioning\n")
 
 F = transpose(T2)*F2*T1
 
endfunction


function F = force_rank2(F) % Force singularity constraint det(F)=0
% ==============
[U, D, V] = svd(F); % Singular value decomposition
D(3, 3) = 0; % Smallest singular value must be 0
F = U * D * V';
endfunction

 
 function [e,ePrime]=calculate_epipoles(F)
   [U,D,V] = svd(F)
   
   e = [V(:,end)]
   ePrime = [U(:,end)]
   
 endfunction

 
 function [Pn,Pprime]=get_projection_matrix(F,ePrime)
    Pn = [1 0 0 0;0 1 0 0;0 0 1 0]
    ePrime2 = [0 -ePrime(3,1) ePrime(2,1); ePrime(3,1) 0 -ePrime(1,1); -ePrime(2,1) ePrime(1,1) 0]
    Pprime1 = [ePrime2*F+[ePrime ePrime ePrime]]
    Pprime = [Pprime1 ePrime]
 endfunction
 
 
  function [xy,xyprime]=take_image_points_from_file()
    fh = fopen('sample.dat', 'r');
    A = fscanf(fh, '%f%f%f%f', [4 inf]);
    fclose(fh);
    xy = cell(1,2)
    xyprime = cell(1,2)
    for k = 1:length(A)
    xy{k} = A(1:2,k)
    xyprime{k} = A(3:4,k)
    end
 

endfunction


function [xy_cp,xyprime_cp,XYZ1,XYZ2,XYZ3,XYZ4,XYZ5]=take_control_points_from_file()
    fh = fopen('pp.dat', 'r');
    A = fscanf(fh, '%f%f%f%f%f%f%f', [7 inf]);
    fclose(fh);
    xy_cp = cell(1,2)
    xyprime_cp = cell(1,2)
    
    for k = 1:5
      xy_cp{k} = A(1:2,k)
      xyprime_cp{k} = A(3:4,k)
    end

    XYZ1 = A(5:7,1)
    XYZ2 = A(5:7,2)
    XYZ3 = A(5:7,3)
    XYZ4 = A(5:7,4)
    XYZ5 = A(5:7,5)
  endfunction
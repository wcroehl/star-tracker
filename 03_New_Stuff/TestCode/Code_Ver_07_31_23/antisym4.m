function A = antisym4(vec) 
A = [      0   vec(3) -vec(2)  vec(1)
      -vec(3)      0   vec(1)  vec(2)   
       vec(2) -vec(1)      0   vec(3)
      -vec(1) -vec(2) -vec(3)      0 ];    
end
function [ Dmax ] = FindW( distgeo, maxdist,Ax, Ay, Az, Shape,R )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%%
%define initial Dmax

Dmax = zeros(Ax,Ay,Az);
w =  zeros(Ax,Ay,Az);
Dmax= double(distgeo);
%Finding Dmax
%iterate rounds
for t = 1: round(maxdist)  
    disp(['GPA recursion reached ',num2str(t/round(maxdist)*100),'%'])
    %for i[nternal cells
    for i = 2:Ax-1
        for j = 2:Ay-1
            for k = 2:Az-1
                %if location is pore
                if distgeo(i,j,k) ~= 0
                    %check surroundings and assign max as 
                    if distgeo(i,j,k) < distgeo(i+1,j,k)
                      if Dmax(i,j,k) < Dmax(i+1,j,k)
                      Dmax(i,j,k) = Dmax(i+1,j,k);
                      end
                    end

                    if distgeo(i,j,k) < distgeo(i-1,j,k)
                      if Dmax(i,j,k) < Dmax(i-1,j,k)
                          Dmax(i,j,k) = Dmax(i-1,j,k);
                      end
                    end 

                    if distgeo(i,j,k) < distgeo(i,j+1,k)
                      if Dmax(i,j,k) < Dmax(i,j+1,k)
                          Dmax(i,j,k) = Dmax(i,j+1,k);
                      end
                    end 

                    if distgeo(i,j,k) < distgeo(i,j-1,k)
                      if Dmax(i,j,k) < Dmax(i,j-1,k)
                          Dmax(i,j,k) = Dmax(i,j-1,k);
                      end
                    end

                    if distgeo(i,j,k) < distgeo(i,j,k+1)
                      if Dmax(i,j,k) < Dmax(i,j,k+1)
                          Dmax(i,j,k) = Dmax(i,j,k+1);
                      end
                    end

                    if distgeo(i,j,k) < distgeo(i,j,k-1)
                      if Dmax(i,j,k) < Dmax(i,j,k-1)
                         Dmax(i,j,k) = Dmax(i,j,k-1);
                      end
                    end


                    if distgeo(i,j,k) < distgeo(i-1,j+1,k+1)
                      if Dmax(i,j,k) < Dmax(i-1,j+1,k+1)
                          Dmax(i,j,k) = Dmax(i-1,j+1,k+1);
                      end
                    end

                     if distgeo(i,j,k) < distgeo(i,j+1,k+1)
                      if Dmax(i,j,k) < Dmax(i,j+1,k+1)
                          Dmax(i,j,k) = Dmax(i,j+1,k+1);
                      end
                    end

                    if distgeo(i,j,k) < distgeo(i+1,j+1,k+1)
                      if Dmax(i,j,k) < Dmax(i+1,j+1,k+1)
                          Dmax(i,j,k) = Dmax(i+1,j+1,k+1);
                      end
                    end

                    if distgeo(i,j,k) < distgeo(i-1,j,k+1)
                      if Dmax(i,j,k) < Dmax(i-1,j,k+1)
                          Dmax(i,j,k) = Dmax(i-1,j,k+1);
                      end
                    end

                    if distgeo(i,j,k) < distgeo(i+1,j,k+1)
                      if Dmax(i,j,k) < Dmax(i+1,j,k+1)
                          Dmax(i,j,k) = Dmax(i+1,j,k+1);
                      end
                    end

                    if distgeo(i,j,k) < distgeo(i-1,j-1,k+1)
                      if Dmax(i,j,k) < Dmax(i-1,j-1,k+1)
                          Dmax(i,j,k) = Dmax(i-1,j-1,k+1);
                      end
                    end

                    if distgeo(i,j,k) < distgeo(i,j-1,k+1)
                      if Dmax(i,j,k) < Dmax(i,j-1,k+1)
                          Dmax(i,j,k) = Dmax(i,j-1,k+1);
                      end
                    end

                    if distgeo(i,j,k) < distgeo(i+1,j-1,k+1)
                      if Dmax(i,j,k) < Dmax(i+1,j-1,k+1)
                          Dmax(i,j,k) = Dmax(i+1,j-1,k+1);
                      end
                    end

                    if distgeo(i,j,k) < distgeo(i-1,j+1,k)
                      if Dmax(i,j,k) < Dmax(i-1,j+1,k)
                      Dmax(i,j,k) = Dmax(i-1,j+1,k);
                      end
                    end

                    if distgeo(i,j,k) < distgeo(i+1,j+1,k)
                      if Dmax(i,j,k) < Dmax(i+1,j+1,k)
                      Dmax(i,j,k) = Dmax(i+1,j+1,k);
                      end
                    end

                    if distgeo(i,j,k) < distgeo(i+1,j-1,k)
                      if Dmax(i,j,k) < Dmax(i+1,j-1,k)
                      Dmax(i,j,k) = Dmax(i+1,j-1,k);
                      end
                    end

                    if distgeo(i,j,k) < distgeo(i-1,j-1,k)
                      if Dmax(i,j,k) < Dmax(i-1,j-1,k)
                      Dmax(i,j,k) = Dmax(i-1,j-1,k);
                      end
                    end

                    if distgeo(i,j,k) < distgeo(i-1,j+1,k-1)
                      if Dmax(i,j,k) < Dmax(i-1,j+1,k-1)
                          Dmax(i,j,k) = Dmax(i-1,j+1,k-1);
                      end
                    end

                     if distgeo(i,j,k) < distgeo(i,j+1,k-1)
                      if Dmax(i,j,k) < Dmax(i,j+1,k-1)
                          Dmax(i,j,k) = Dmax(i,j+1,k-1);
                      end
                      end

                    if distgeo(i,j,k) < distgeo(i+1,j+1,k-1)
                      if Dmax(i,j,k) < Dmax(i+1,j+1,k-1)
                          Dmax(i,j,k) = Dmax(i+1,j+1,k-1);
                      end
                    end

                    if distgeo(i,j,k) < distgeo(i-1,j,k-1)
                      if Dmax(i,j,k) < Dmax(i-1,j,k-1)
                          Dmax(i,j,k) = Dmax(i-1,j,k-1);
                      end
                    end

                    if distgeo(i,j,k) < distgeo(i+1,j,k-1)
                      if Dmax(i,j,k) < Dmax(i+1,j,k-1)
                          Dmax(i,j,k) = Dmax(i+1,j,k-1);
                      end
                    end

                    if distgeo(i,j,k) < distgeo(i-1,j-1,k-1)
                      if Dmax(i,j,k) < Dmax(i-1,j-1,k-1)
                          Dmax(i,j,k) = Dmax(i-1,j-1,k-1);
                      end
                    end

                    if distgeo(i,j,k) < distgeo(i,j-1,k-1)
                      if Dmax(i,j,k) < Dmax(i,j-1,k-1)
                          Dmax(i,j,k) = Dmax(i,j-1,k-1);
                      end
                    end

                    if distgeo(i,j,k) < distgeo(i+1,j-1,k-1)
                      if Dmax(i,j,k) < Dmax(i+1,j-1,k-1)
                          Dmax(i,j,k) = Dmax(i+1,j-1,k-1);
                      end
                    end
      
                    end
                   end
                  end
        end
   
%face k = 1

for i = 2:Ax-1
         for j = 2:Ay-1
             
            k = 1;
            
      if distgeo(i,j,k) ~= 0
          
      if distgeo(i,j,k) < distgeo(i+1,j,k)
          if Dmax(i,j,k) < Dmax(i+1,j,k)
          Dmax(i,j,k) = Dmax(i+1,j,k);
          end
      end
      
      if distgeo(i,j,k) < distgeo(i-1,j,k)
          if Dmax(i,j,k) < Dmax(i-1,j,k)
              Dmax(i,j,k) = Dmax(i-1,j,k);
          end
      end 
          
      if distgeo(i,j,k) < distgeo(i,j+1,k)
          if Dmax(i,j,k) < Dmax(i,j+1,k)
              Dmax(i,j,k) = Dmax(i,j+1,k);
          end
      end 
            
      if distgeo(i,j,k) < distgeo(i,j-1,k)
          if Dmax(i,j,k) < Dmax(i,j-1,k)
              Dmax(i,j,k) = Dmax(i,j-1,k);
          end
      end
              
      if distgeo(i,j,k) < distgeo(i,j,k+1)
          if Dmax(i,j,k) < Dmax(i,j,k+1)
              Dmax(i,j,k) = Dmax(i,j,k+1);
          end
      end
              
     
      
      
      if distgeo(i,j,k) < distgeo(i-1,j+1,k+1)
          if Dmax(i,j,k) < Dmax(i-1,j+1,k+1)
              Dmax(i,j,k) = Dmax(i-1,j+1,k+1);
          end
      end
      
         if distgeo(i,j,k) < distgeo(i,j+1,k+1)
          if Dmax(i,j,k) < Dmax(i,j+1,k+1)
              Dmax(i,j,k) = Dmax(i,j+1,k+1);
          end
      end
      
       if distgeo(i,j,k) < distgeo(i+1,j+1,k+1)
          if Dmax(i,j,k) < Dmax(i+1,j+1,k+1)
              Dmax(i,j,k) = Dmax(i+1,j+1,k+1);
          end
       end
      
       if distgeo(i,j,k) < distgeo(i-1,j,k+1)
          if Dmax(i,j,k) < Dmax(i-1,j,k+1)
              Dmax(i,j,k) = Dmax(i-1,j,k+1);
          end
       end
      
       if distgeo(i,j,k) < distgeo(i+1,j,k+1)
          if Dmax(i,j,k) < Dmax(i+1,j,k+1)
              Dmax(i,j,k) = Dmax(i+1,j,k+1);
          end
       end
      
        if distgeo(i,j,k) < distgeo(i-1,j-1,k+1)
          if Dmax(i,j,k) < Dmax(i-1,j-1,k+1)
              Dmax(i,j,k) = Dmax(i-1,j-1,k+1);
          end
        end
       
        if distgeo(i,j,k) < distgeo(i,j-1,k+1)
          if Dmax(i,j,k) < Dmax(i,j-1,k+1)
              Dmax(i,j,k) = Dmax(i,j-1,k+1);
          end
        end
       
        if distgeo(i,j,k) < distgeo(i+1,j-1,k+1)
          if Dmax(i,j,k) < Dmax(i+1,j-1,k+1)
              Dmax(i,j,k) = Dmax(i+1,j-1,k+1);
          end
        end
       
        if distgeo(i,j,k) < distgeo(i-1,j+1,k)
          if Dmax(i,j,k) < Dmax(i-1,j+1,k)
          Dmax(i,j,k) = Dmax(i-1,j+1,k);
          end
        end
      
      if distgeo(i,j,k) < distgeo(i+1,j+1,k)
          if Dmax(i,j,k) < Dmax(i+1,j+1,k)
          Dmax(i,j,k) = Dmax(i+1,j+1,k);
          end
      end
      
      if distgeo(i,j,k) < distgeo(i+1,j-1,k)
          if Dmax(i,j,k) < Dmax(i+1,j-1,k)
          Dmax(i,j,k) = Dmax(i+1,j-1,k);
          end
      end
      
      if distgeo(i,j,k) < distgeo(i-1,j-1,k)
          if Dmax(i,j,k) < Dmax(i-1,j-1,k)
          Dmax(i,j,k) = Dmax(i-1,j-1,k);
          end
      end
  
                    end
                   end
end
                  

%face k = end

for i = 2:Ax-1
         for j = 2:Ay-1
             
            k = Az;
            
      if distgeo(i,j,k) ~= 0
          
      if distgeo(i,j,k) < distgeo(i+1,j,k)
          if Dmax(i,j,k) < Dmax(i+1,j,k)
          Dmax(i,j,k) = Dmax(i+1,j,k);
          end
      end
      
      if distgeo(i,j,k) < distgeo(i-1,j,k)
          if Dmax(i,j,k) < Dmax(i-1,j,k)
              Dmax(i,j,k) = Dmax(i-1,j,k);
          end
      end 
          
      if distgeo(i,j,k) < distgeo(i,j+1,k)
          if Dmax(i,j,k) < Dmax(i,j+1,k)
              Dmax(i,j,k) = Dmax(i,j+1,k);
          end
      end 
            
      if distgeo(i,j,k) < distgeo(i,j-1,k)
          if Dmax(i,j,k) < Dmax(i,j-1,k)
              Dmax(i,j,k) = Dmax(i,j-1,k);
          end
      end
              
      if distgeo(i,j,k) < distgeo(i,j,k-1)
          if Dmax(i,j,k) < Dmax(i,j,k-1)
              Dmax(i,j,k) = Dmax(i,j,k-1);
          end
      end
              
     
      
      
      if distgeo(i,j,k) < distgeo(i-1,j+1,k-1)
          if Dmax(i,j,k) < Dmax(i-1,j+1,k-1)
              Dmax(i,j,k) = Dmax(i-1,j+1,k-1);
          end
      end
      
         if distgeo(i,j,k) < distgeo(i,j+1,k-1)
          if Dmax(i,j,k) < Dmax(i,j+1,k-1)
              Dmax(i,j,k) = Dmax(i,j+1,k-1);
          end
      end
      
       if distgeo(i,j,k) < distgeo(i+1,j+1,k-1)
          if Dmax(i,j,k) < Dmax(i+1,j+1,k-1)
              Dmax(i,j,k) = Dmax(i+1,j+1,k-1);
          end
       end
      
       if distgeo(i,j,k) < distgeo(i-1,j,k-1)
          if Dmax(i,j,k) < Dmax(i-1,j,k-1)
              Dmax(i,j,k) = Dmax(i-1,j,k-1);
          end
       end
      
       if distgeo(i,j,k) < distgeo(i+1,j,k-1)
          if Dmax(i,j,k) < Dmax(i+1,j,k-1)
              Dmax(i,j,k) = Dmax(i+1,j,k-1);
          end
       end
      
        if distgeo(i,j,k) < distgeo(i-1,j-1,k-1)
          if Dmax(i,j,k) < Dmax(i-1,j-1,k-1)
              Dmax(i,j,k) = Dmax(i-1,j-1,k-1);
          end
        end
       
        if distgeo(i,j,k) < distgeo(i,j-1,k-1)
          if Dmax(i,j,k) < Dmax(i,j-1,k-1)
              Dmax(i,j,k) = Dmax(i,j-1,k-1);
          end
        end
       
        if distgeo(i,j,k) < distgeo(i+1,j-1,k-1)
          if Dmax(i,j,k) < Dmax(i+1,j-1,k-1)
              Dmax(i,j,k) = Dmax(i+1,j-1,k-1);
          end
        end
       
        if distgeo(i,j,k) < distgeo(i-1,j+1,k)
          if Dmax(i,j,k) < Dmax(i-1,j+1,k)
          Dmax(i,j,k) = Dmax(i-1,j+1,k);
          end
        end
      
      if distgeo(i,j,k) < distgeo(i+1,j+1,k)
          if Dmax(i,j,k) < Dmax(i+1,j+1,k)
          Dmax(i,j,k) = Dmax(i+1,j+1,k);
          end
      end
      
      if distgeo(i,j,k) < distgeo(i+1,j-1,k)
          if Dmax(i,j,k) < Dmax(i+1,j-1,k)
          Dmax(i,j,k) = Dmax(i+1,j-1,k);
          end
      end
      
      if distgeo(i,j,k) < distgeo(i-1,j-1,k)
          if Dmax(i,j,k) < Dmax(i-1,j-1,k)
          Dmax(i,j,k) = Dmax(i-1,j-1,k);
          end
      end
  
                    end
                   end
end
       
                  
%face j = 1

  for i = 2:Ax-1
         
            for k = 2:Az-1
                
                j = 1;
            
      if distgeo(i,j,k) ~= 0
          
      if distgeo(i,j,k) < distgeo(i+1,j,k)
          if Dmax(i,j,k) < Dmax(i+1,j,k)
          Dmax(i,j,k) = Dmax(i+1,j,k);
          end
      end
      
      if distgeo(i,j,k) < distgeo(i-1,j,k)
          if Dmax(i,j,k) < Dmax(i-1,j,k)
              Dmax(i,j,k) = Dmax(i-1,j,k);
          end
      end 
          
      if distgeo(i,j,k) < distgeo(i,j+1,k)
          if Dmax(i,j,k) < Dmax(i,j+1,k)
              Dmax(i,j,k) = Dmax(i,j+1,k);
          end
      end 
            
     
              
      if distgeo(i,j,k) < distgeo(i,j,k+1)
          if Dmax(i,j,k) < Dmax(i,j,k+1)
              Dmax(i,j,k) = Dmax(i,j,k+1);
          end
      end
              
      if distgeo(i,j,k) < distgeo(i,j,k-1)
          if Dmax(i,j,k) < Dmax(i,j,k-1)
             Dmax(i,j,k) = Dmax(i,j,k-1);
          end
      end
      
      
      if distgeo(i,j,k) < distgeo(i-1,j+1,k+1)
          if Dmax(i,j,k) < Dmax(i-1,j+1,k+1)
              Dmax(i,j,k) = Dmax(i-1,j+1,k+1);
          end
      end
      
         if distgeo(i,j,k) < distgeo(i,j+1,k+1)
          if Dmax(i,j,k) < Dmax(i,j+1,k+1)
              Dmax(i,j,k) = Dmax(i,j+1,k+1);
          end
      end
      
       if distgeo(i,j,k) < distgeo(i+1,j+1,k+1)
          if Dmax(i,j,k) < Dmax(i+1,j+1,k+1)
              Dmax(i,j,k) = Dmax(i+1,j+1,k+1);
          end
       end
      
       if distgeo(i,j,k) < distgeo(i-1,j,k+1)
          if Dmax(i,j,k) < Dmax(i-1,j,k+1)
              Dmax(i,j,k) = Dmax(i-1,j,k+1);
          end
       end
      
       if distgeo(i,j,k) < distgeo(i+1,j,k+1)
          if Dmax(i,j,k) < Dmax(i+1,j,k+1)
              Dmax(i,j,k) = Dmax(i+1,j,k+1);
          end
       end
      
        
       
        if distgeo(i,j,k) < distgeo(i-1,j+1,k)
          if Dmax(i,j,k) < Dmax(i-1,j+1,k)
          Dmax(i,j,k) = Dmax(i-1,j+1,k);
          end
        end
      
      if distgeo(i,j,k) < distgeo(i+1,j+1,k)
          if Dmax(i,j,k) < Dmax(i+1,j+1,k)
          Dmax(i,j,k) = Dmax(i+1,j+1,k);
          end
      end
      
     
      
      if distgeo(i,j,k) < distgeo(i-1,j+1,k-1)
          if Dmax(i,j,k) < Dmax(i-1,j+1,k-1)
              Dmax(i,j,k) = Dmax(i-1,j+1,k-1);
          end
      end
      
         if distgeo(i,j,k) < distgeo(i,j+1,k-1)
          if Dmax(i,j,k) < Dmax(i,j+1,k-1)
              Dmax(i,j,k) = Dmax(i,j+1,k-1);
          end
          end
      
       if distgeo(i,j,k) < distgeo(i+1,j+1,k-1)
          if Dmax(i,j,k) < Dmax(i+1,j+1,k-1)
              Dmax(i,j,k) = Dmax(i+1,j+1,k-1);
          end
       end
      
       if distgeo(i,j,k) < distgeo(i-1,j,k-1)
          if Dmax(i,j,k) < Dmax(i-1,j,k-1)
              Dmax(i,j,k) = Dmax(i-1,j,k-1);
          end
       end
      
       if distgeo(i,j,k) < distgeo(i+1,j,k-1)
          if Dmax(i,j,k) < Dmax(i+1,j,k-1)
              Dmax(i,j,k) = Dmax(i+1,j,k-1);
          end
       end
      
                    end
                   end
                  end
      
%face j = end

  for i = 2:Ax-1
       
  for k = 2:Az-1
                
              j = Ay;
      if distgeo(i,j,k) ~= 0
          
      if distgeo(i,j,k) < distgeo(i+1,j,k)
          if Dmax(i,j,k) < Dmax(i+1,j,k)
          Dmax(i,j,k) = Dmax(i+1,j,k);
          end
      end
      
      if distgeo(i,j,k) < distgeo(i-1,j,k)
          if Dmax(i,j,k) < Dmax(i-1,j,k)
              Dmax(i,j,k) = Dmax(i-1,j,k);
          end
      end 
          
      if distgeo(i,j,k) < distgeo(i,j-1,k)
          if Dmax(i,j,k) < Dmax(i,j-1,k)
              Dmax(i,j,k) = Dmax(i,j-1,k);
          end
      end 
            
     
              
      if distgeo(i,j,k) < distgeo(i,j,k+1)
          if Dmax(i,j,k) < Dmax(i,j,k+1)
              Dmax(i,j,k) = Dmax(i,j,k+1);
          end
      end
              
      if distgeo(i,j,k) < distgeo(i,j,k-1)
          if Dmax(i,j,k) < Dmax(i,j,k-1)
             Dmax(i,j,k) = Dmax(i,j,k-1);
          end
      end
      
      
      if distgeo(i,j,k) < distgeo(i-1,j-1,k+1)
          if Dmax(i,j,k) < Dmax(i-1,j-1,k+1)
              Dmax(i,j,k) = Dmax(i-1,j-1,k+1);
          end
      end
      
         if distgeo(i,j,k) < distgeo(i,j-1,k+1)
          if Dmax(i,j,k) < Dmax(i,j-1,k+1)
              Dmax(i,j,k) = Dmax(i,j-1,k+1);
          end
      end
      
       if distgeo(i,j,k) < distgeo(i+1,j-1,k+1)
          if Dmax(i,j,k) < Dmax(i+1,j-1,k+1)
              Dmax(i,j,k) = Dmax(i+1,j-1,k+1);
          end
       end
      
       if distgeo(i,j,k) < distgeo(i-1,j,k+1)
          if Dmax(i,j,k) < Dmax(i-1,j,k+1)
              Dmax(i,j,k) = Dmax(i-1,j,k+1);
          end
       end
      
       if distgeo(i,j,k) < distgeo(i+1,j,k+1)
          if Dmax(i,j,k) < Dmax(i+1,j,k+1)
              Dmax(i,j,k) = Dmax(i+1,j,k+1);
          end
       end
      
        
       
        if distgeo(i,j,k) < distgeo(i-1,j-1,k)
          if Dmax(i,j,k) < Dmax(i-1,j-1,k)
          Dmax(i,j,k) = Dmax(i-1,j-1,k);
          end
        end
      
      if distgeo(i,j,k) < distgeo(i+1,j-1,k)
          if Dmax(i,j,k) < Dmax(i+1,j-1,k)
          Dmax(i,j,k) = Dmax(i+1,j-1,k);
          end
      end
      
     
      
      if distgeo(i,j,k) < distgeo(i-1,j-1,k-1)
          if Dmax(i,j,k) < Dmax(i-1,j-1,k-1)
              Dmax(i,j,k) = Dmax(i-1,j-1,k-1);
          end
      end
      
         if distgeo(i,j,k) < distgeo(i,j-1,k-1)
          if Dmax(i,j,k) < Dmax(i,j-1,k-1)
              Dmax(i,j,k) = Dmax(i,j-1,k-1);
          end
          end
      
       if distgeo(i,j,k) < distgeo(i+1,j-1,k-1)
          if Dmax(i,j,k) < Dmax(i+1,j-1,k-1)
              Dmax(i,j,k) = Dmax(i+1,j-1,k-1);
          end
       end
      
       if distgeo(i,j,k) < distgeo(i-1,j,k-1)
          if Dmax(i,j,k) < Dmax(i-1,j,k-1)
              Dmax(i,j,k) = Dmax(i-1,j,k-1);
          end
       end
      
       if distgeo(i,j,k) < distgeo(i+1,j,k-1)
          if Dmax(i,j,k) < Dmax(i+1,j,k-1)
              Dmax(i,j,k) = Dmax(i+1,j,k-1);
          end
       end
      
                    end
                   end
  end
                  
  
%face i = 1

         for j = 2:Ay-1
            for k = 2:Az-1
                
                i = 1;
            
      if distgeo(i,j,k) ~= 0
          
      if distgeo(i,j,k) < distgeo(i+1,j,k)
          if Dmax(i,j,k) < Dmax(i+1,j,k)
          Dmax(i,j,k) = Dmax(i+1,j,k);
          end
      end
      
     
          
      if distgeo(i,j,k) < distgeo(i,j+1,k)
          if Dmax(i,j,k) < Dmax(i,j+1,k)
              Dmax(i,j,k) = Dmax(i,j+1,k);
          end
      end 
            
      if distgeo(i,j,k) < distgeo(i,j-1,k)
          if Dmax(i,j,k) < Dmax(i,j-1,k)
              Dmax(i,j,k) = Dmax(i,j-1,k);
          end
      end
              
      if distgeo(i,j,k) < distgeo(i,j,k+1)
          if Dmax(i,j,k) < Dmax(i,j,k+1)
              Dmax(i,j,k) = Dmax(i,j,k+1);
          end
      end
              
      if distgeo(i,j,k) < distgeo(i,j,k-1)
          if Dmax(i,j,k) < Dmax(i,j,k-1)
             Dmax(i,j,k) = Dmax(i,j,k-1);
          end
      end
      
      
      
      
         if distgeo(i,j,k) < distgeo(i,j+1,k+1)
          if Dmax(i,j,k) < Dmax(i,j+1,k+1)
              Dmax(i,j,k) = Dmax(i,j+1,k+1);
          end
      end
      
       if distgeo(i,j,k) < distgeo(i+1,j+1,k+1)
          if Dmax(i,j,k) < Dmax(i+1,j+1,k+1)
              Dmax(i,j,k) = Dmax(i+1,j+1,k+1);
          end
       end
      
    
      
       if distgeo(i,j,k) < distgeo(i+1,j,k+1)
          if Dmax(i,j,k) < Dmax(i+1,j,k+1)
              Dmax(i,j,k) = Dmax(i+1,j,k+1);
          end
       end
      
       
       
        if distgeo(i,j,k) < distgeo(i,j-1,k+1)
          if Dmax(i,j,k) < Dmax(i,j-1,k+1)
              Dmax(i,j,k) = Dmax(i,j-1,k+1);
          end
        end
       
        if distgeo(i,j,k) < distgeo(i+1,j-1,k+1)
          if Dmax(i,j,k) < Dmax(i+1,j-1,k+1)
              Dmax(i,j,k) = Dmax(i+1,j-1,k+1);
          end
        end
       
     
      if distgeo(i,j,k) < distgeo(i+1,j+1,k)
          if Dmax(i,j,k) < Dmax(i+1,j+1,k)
          Dmax(i,j,k) = Dmax(i+1,j+1,k);
          end
      end
      
      if distgeo(i,j,k) < distgeo(i+1,j-1,k)
          if Dmax(i,j,k) < Dmax(i+1,j-1,k)
          Dmax(i,j,k) = Dmax(i+1,j-1,k);
          end
      end
      
      
      
         if distgeo(i,j,k) < distgeo(i,j+1,k-1)
          if Dmax(i,j,k) < Dmax(i,j+1,k-1)
              Dmax(i,j,k) = Dmax(i,j+1,k-1);
          end
          end
      
       if distgeo(i,j,k) < distgeo(i+1,j+1,k-1)
          if Dmax(i,j,k) < Dmax(i+1,j+1,k-1)
              Dmax(i,j,k) = Dmax(i+1,j+1,k-1);
          end
       end
      
     
      
       if distgeo(i,j,k) < distgeo(i+1,j,k-1)
          if Dmax(i,j,k) < Dmax(i+1,j,k-1)
              Dmax(i,j,k) = Dmax(i+1,j,k-1);
          end
       end
      
       
       
        if distgeo(i,j,k) < distgeo(i,j-1,k-1)
          if Dmax(i,j,k) < Dmax(i,j-1,k-1)
              Dmax(i,j,k) = Dmax(i,j-1,k-1);
          end
        end
       
        if distgeo(i,j,k) < distgeo(i+1,j-1,k-1)
          if Dmax(i,j,k) < Dmax(i+1,j-1,k-1)
              Dmax(i,j,k) = Dmax(i+1,j-1,k-1);
          end
        end
      
                    end
                   end
         end
                  
         
%face i = end

         for j = 2:Ay-1
            for k = 2:Az-1
                
                i = Ax;
            
      if distgeo(i,j,k) ~= 0
          
      if distgeo(i,j,k) < distgeo(i-1,j,k)
          if Dmax(i,j,k) < Dmax(i-1,j,k)
          Dmax(i,j,k) = Dmax(i-1,j,k);
          end
      end
      
     
          
      if distgeo(i,j,k) < distgeo(i,j+1,k)
          if Dmax(i,j,k) < Dmax(i,j+1,k)
              Dmax(i,j,k) = Dmax(i,j+1,k);
          end
      end 
            
      if distgeo(i,j,k) < distgeo(i,j-1,k)
          if Dmax(i,j,k) < Dmax(i,j-1,k)
              Dmax(i,j,k) = Dmax(i,j-1,k);
          end
      end
              
      if distgeo(i,j,k) < distgeo(i,j,k+1)
          if Dmax(i,j,k) < Dmax(i,j,k+1)
              Dmax(i,j,k) = Dmax(i,j,k+1);
          end
      end
              
      if distgeo(i,j,k) < distgeo(i,j,k-1)
          if Dmax(i,j,k) < Dmax(i,j,k-1)
             Dmax(i,j,k) = Dmax(i,j,k-1);
          end
      end
      
      
      
      
         if distgeo(i,j,k) < distgeo(i,j+1,k+1)
          if Dmax(i,j,k) < Dmax(i,j+1,k+1)
              Dmax(i,j,k) = Dmax(i,j+1,k+1);
          end
      end
      
       if distgeo(i,j,k) < distgeo(i-1,j+1,k+1)
          if Dmax(i,j,k) < Dmax(i-1,j+1,k+1)
              Dmax(i,j,k) = Dmax(i-1,j+1,k+1);
          end
       end
      
    
      
       if distgeo(i,j,k) < distgeo(i-1,j,k+1)
          if Dmax(i,j,k) < Dmax(i-1,j,k+1)
              Dmax(i,j,k) = Dmax(i-1,j,k+1);
          end
       end
      
       
       
        if distgeo(i,j,k) < distgeo(i,j-1,k+1)
          if Dmax(i,j,k) < Dmax(i,j-1,k+1)
              Dmax(i,j,k) = Dmax(i,j-1,k+1);
          end
        end
       
        if distgeo(i,j,k) < distgeo(i-1,j-1,k+1)
          if Dmax(i,j,k) < Dmax(i-1,j-1,k+1)
              Dmax(i,j,k) = Dmax(i-1,j-1,k+1);
          end
        end
       
     
      if distgeo(i,j,k) < distgeo(i-1,j+1,k)
          if Dmax(i,j,k) < Dmax(i-1,j+1,k)
          Dmax(i,j,k) = Dmax(i-1,j+1,k);
          end
      end
      
      if distgeo(i,j,k) < distgeo(i-1,j-1,k)
          if Dmax(i,j,k) < Dmax(i-1,j-1,k)
          Dmax(i,j,k) = Dmax(i-1,j-1,k);
          end
      end
      
      
      
         if distgeo(i,j,k) < distgeo(i,j+1,k-1)
          if Dmax(i,j,k) < Dmax(i,j+1,k-1)
              Dmax(i,j,k) = Dmax(i,j+1,k-1);
          end
          end
      
       if distgeo(i,j,k) < distgeo(i-1,j+1,k-1)
          if Dmax(i,j,k) < Dmax(i-1,j+1,k-1)
              Dmax(i,j,k) = Dmax(i-1,j+1,k-1);
          end
       end
      
     
      
       if distgeo(i,j,k) < distgeo(i-1,j,k-1)
          if Dmax(i,j,k) < Dmax(i-1,j,k-1)
              Dmax(i,j,k) = Dmax(i-1,j,k-1);
          end
       end
      
       
       
        if distgeo(i,j,k) < distgeo(i,j-1,k-1)
          if Dmax(i,j,k) < Dmax(i,j-1,k-1)
              Dmax(i,j,k) = Dmax(i,j-1,k-1);
          end
        end
       
        if distgeo(i,j,k) < distgeo(i-1,j-1,k-1)
          if Dmax(i,j,k) < Dmax(i-1,j-1,k-1)
              Dmax(i,j,k) = Dmax(i-1,j-1,k-1);
          end
        end
      
                    end
                   end
                  end
       
       
  
  
        
  %check i edge
  for i=2:Ax-1
      
      if distgeo(i,j,k) ~= 0
  
    if distgeo(i,1,1) < distgeo(i+1,1,1)
       if Dmax(i,1,1) < Dmax(i+1,1,1)
            Dmax(i,1,1) = Dmax(i+1,1,1);
       end
   end
   
   if distgeo(i,1,1) < distgeo(i-1,1,1)
       if Dmax(i,1,1) < Dmax(i-1,1,1)
            Dmax(i,1,1) = Dmax(i-1,1,1);
       end
   end
   
    if distgeo(i,1,1) < distgeo(i-1,2,1)
       if Dmax(i,1,1) < Dmax(i-1,2,1)
            Dmax(i,1,1) = Dmax(i-1,2,1);
       end
    end
   
     if distgeo(i,1,1) < distgeo(i+1,2,1)
       if Dmax(i,1,1) < Dmax(i+1,2,1)
            Dmax(i,1,1) = Dmax(i+1,2,1);
       end
     end
   
     if distgeo(i,1,1) < distgeo(i,2,1)
       if Dmax(i,1,1) < Dmax(i,2,1)
            Dmax(i,1,1) = Dmax(i,2,1);
       end
     end
   
       if distgeo(i,1,1) < distgeo(i,1,2)
       if Dmax(i,1,1) < Dmax(i,1,2)
            Dmax(i,1,1) = Dmax(i,1,2);
       end
      end
      if distgeo(i,1,1) < distgeo(i+1,1,2)
       if Dmax(i,1,1) < Dmax(i+1,1,2)
            Dmax(i,1,1) = Dmax(i+1,1,2);
       end
   end
   
   if distgeo(i,1,1) < distgeo(i-1,1,2)
       if Dmax(i,1,1) < Dmax(i-1,1,2)
            Dmax(i,1,1) = Dmax(i-1,1,2);
       end
   end
   
    if distgeo(i,1,1) < distgeo(i-1,2,2)
       if Dmax(i,1,1) < Dmax(i-1,2,2)
            Dmax(i,1,1) = Dmax(i-1,2,2);
       end
    end
   
     if distgeo(i,1,1) < distgeo(i+1,2,2)
       if Dmax(i,1,1) < Dmax(i+1,2,2)
            Dmax(i,1,1) = Dmax(i+1,2,2);
       end
     end
   
     if distgeo(i,1,1) < distgeo(i,2,2)
       if Dmax(i,1,1) < Dmax(i,2,2)
            Dmax(i,1,1) = Dmax(i,2,2);
       end
     end
   
 %%
 
     if distgeo(i,end,1) < distgeo(i+1,end,1)
       if Dmax(i,end,1) < Dmax(i+1,end,1)
            Dmax(i,end,1) = Dmax(i+1,end,1);
       end
   end
   
   if distgeo(i,end,1) < distgeo(i-1,end,1)
       if Dmax(i,end,1) < Dmax(i-1,end,1)
            Dmax(i,end,1) = Dmax(i-1,end,1);
       end
   end
   
    if distgeo(i,end,1) < distgeo(i-1,end-1,1)
       if Dmax(i,end,1) < Dmax(i-1,end-1,1)
            Dmax(i,end,1) = Dmax(i-1,end-1,1);
       end
    end
   
     if distgeo(i,end,1) < distgeo(i+1,end-1,1)
       if Dmax(i,end,1) < Dmax(i+1,end-1,1)
            Dmax(i,end,1) = Dmax(i+1,end-1,1);
       end
     end
   
     if distgeo(i,end,1) < distgeo(i,end-1,1)
       if Dmax(i,end,1) < Dmax(i,end-1,1)
            Dmax(i,end,1) = Dmax(i,end-1,1);
       end
     end
   
      if distgeo(i,end,1) < distgeo(i,end,2)
       if Dmax(i,end,1) < Dmax(i,end,2)
            Dmax(i,end,1) = Dmax(i,end,2);
       end
      end
   
      if distgeo(i,end,1) < distgeo(i+1,end,2)
       if Dmax(i,end,1) < Dmax(i+1,end,2)
            Dmax(i,end,1) = Dmax(i+1,end,2);
       end
   end
   
   if distgeo(i,end,1) < distgeo(i-1,end,2)
       if Dmax(i,end,1) < Dmax(i-1,end,2)
            Dmax(i,end,1) = Dmax(i-1,end,2);
       end
   end
   
    if distgeo(i,end,1) < distgeo(i-1,end-1,2)
       if Dmax(i,end,1) < Dmax(i-1,end-1,2)
            Dmax(i,end,1) = Dmax(i-1,end-1,2);
       end
    end
   
     if distgeo(i,end,1) < distgeo(i+1,end-1,2)
       if Dmax(i,end,1) < Dmax(i+1,end-1,2)
            Dmax(i,end,1) = Dmax(i+1,end-1,2);
       end
     end
   
     if distgeo(i,end,1) < distgeo(i,end-1,2)
       if Dmax(i,end,1) < Dmax(i,end-1,2)
            Dmax(i,end,1) = Dmax(i,end-1,2);
       end
     end
   
     %% change in z
     
     if distgeo(i,1,end) < distgeo(i+1,1,end)
       if Dmax(i,1,end) < Dmax(i+1,1,end)
            Dmax(i,1,end) = Dmax(i+1,1,end);
       end
   end
   
   if distgeo(i,1,end) < distgeo(i-1,1,end)
       if Dmax(i,1,end) < Dmax(i-1,1,end)
            Dmax(i,1,end) = Dmax(i-1,1,end);
       end
   end
   
    if distgeo(i,1,end) < distgeo(i-1,2,end)
       if Dmax(i,1,end) < Dmax(i-1,2,end)
            Dmax(i,1,end) = Dmax(i-1,2,end);
       end
    end
   
     if distgeo(i,1,end) < distgeo(i+1,2,end)
       if Dmax(i,1,end) < Dmax(i+1,2,end)
            Dmax(i,1,end) = Dmax(i+1,2,end);
       end
     end
   
     if distgeo(i,1,end) < distgeo(i,2,end)
       if Dmax(i,1,end) < Dmax(i,2,end)
            Dmax(i,1,end) = Dmax(i,2,end);
       end
     end
   
       if distgeo(i,end,end) < distgeo(i,1,end-1)
       if Dmax(i,end,end) < Dmax(i,1,end-1)
            Dmax(i,end,end) = Dmax(i,1,end-1);
       end
      end
      if distgeo(i,1,end) < distgeo(i+1,1,end-1)
       if Dmax(i,1,end) < Dmax(i+1,1,end-1)
            Dmax(i,1,end) = Dmax(i+1,1,end-1);
       end
   end
   
   if distgeo(i,1,end) < distgeo(i-1,1,end-1)
       if Dmax(i,1,end) < Dmax(i-1,1,end-1)
            Dmax(i,1,end) = Dmax(i-1,1,end-1);
       end
   end
   
    if distgeo(i,1,end) < distgeo(i-1,2,end-1)
       if Dmax(i,1,end) < Dmax(i-1,2,end-1)
            Dmax(i,1,end) = Dmax(i-1,2,end-1);
       end
    end
   
     if distgeo(i,1,end) < distgeo(i+1,2,end-1)
       if Dmax(i,1,end) < Dmax(i+1,2,end-1)
            Dmax(i,1,end) = Dmax(i+1,2,end-1);
       end
     end
   
     if distgeo(i,1,end) < distgeo(i,2,end-1)
       if Dmax(i,1,end) < Dmax(i,2,end-1)
            Dmax(i,1,end) = Dmax(i,2,end-1);
       end
     end
   
 %%
 
     if distgeo(i,end,end) < distgeo(i+1,end,end)
       if Dmax(i,end,end) < Dmax(i+1,end,end)
            Dmax(i,end,end) = Dmax(i+1,end,end);
       end
   end
   
   if distgeo(i,end,end) < distgeo(i-1,end,end)
       if Dmax(i,end,end) < Dmax(i-1,end,end)
            Dmax(i,end,end) = Dmax(i-1,end,end);
       end
   end
   
    if distgeo(i,end,end) < distgeo(i-1,end-1,end)
       if Dmax(i,end,end) < Dmax(i-1,end-1,end)
            Dmax(i,end,end) = Dmax(i-1,end-1,end);
       end
    end
   
     if distgeo(i,end,end) < distgeo(i+1,end-1,end)
       if Dmax(i,end,end) < Dmax(i+1,end-1,end)
            Dmax(i,end,end) = Dmax(i+1,end-1,end);
       end
     end
   
     if distgeo(i,end,end) < distgeo(i,end-1,end)
       if Dmax(i,end,end) < Dmax(i,end-1,end)
            Dmax(i,end,end) = Dmax(i,end-1,end);
       end
     end
   
      if distgeo(i,end,end) < distgeo(i,end,end-1)
       if Dmax(i,end,end) < Dmax(i,end,end-1)
            Dmax(i,end,end) = Dmax(i,end,end-1);
       end
      end
   
      if distgeo(i,end,end) < distgeo(i+1,end,end-1)
       if Dmax(i,end,end) < Dmax(i+1,end,end-1)
            Dmax(i,end,end) = Dmax(i+1,end,end-1);
       end
   end
   
   if distgeo(i,end,end) < distgeo(i-1,end,end-1)
       if Dmax(i,end,end) < Dmax(i-1,end,end-1)
            Dmax(i,end,end) = Dmax(i-1,end,end-1);
       end
   end
   
    if distgeo(i,end,end) < distgeo(i-1,end-1,end-1)
       if Dmax(i,end,end) < Dmax(i-1,end-1,end-1)
            Dmax(i,end,end) = Dmax(i-1,end-1,end-1);
       end
    end
   
     if distgeo(i,end,end) < distgeo(i+1,end-1,end-1)
       if Dmax(i,end,end) < Dmax(i+1,end-1,end-1)
            Dmax(i,end,end) = Dmax(i+1,end-1,end-1);
       end
     end
   
     if distgeo(i,end,end) < distgeo(i,end-1,end-1)
       if Dmax(i,end,end) < Dmax(i,end-1,end-1)
            Dmax(i,end,end) = Dmax(i,end-1,end-1);
       end
     end
     
      end
  end
     
  %check j edge
  for j=2:Ay-1
  
       if distgeo(i,j,k) ~= 0
    if distgeo(1,j,1) < distgeo(1,j+1,1)
       if Dmax(1,j,1) < Dmax(1,j+1,1)
            Dmax(1,j,1) = Dmax(1,j+1,1);
       end
   end
   
   if distgeo(1,j,1) < distgeo(1,j-1,1)
       if Dmax(1,j,1) < Dmax(1,j-1,1)
            Dmax(1,j,1) = Dmax(1,j-1,1);
       end
   end
   
    if distgeo(1,j,1) < distgeo(2,j-1,1)
       if Dmax(1,j,1) < Dmax(2,j-1,1)
            Dmax(1,j,1) = Dmax(2,j-1,1);
       end
    end
   
     if distgeo(1,j,1) < distgeo(2,j+1,1)
       if Dmax(1,j,1) < Dmax(2,j+1,1)
            Dmax(1,j,1) = Dmax(2,j+1,1);
       end
     end
   
     if distgeo(1,j,1) < distgeo(2,j,1)
       if Dmax(1,j,1) < Dmax(2,j,1)
            Dmax(1,j,1) = Dmax(2,j,1);
       end
     end
   
       if distgeo(1,j,1) < distgeo(1,j,2)
       if Dmax(1,j,1) < Dmax(1,j,2)
            Dmax(1,j,1) = Dmax(1,j,2);
       end
      end
      if distgeo(1,j,1) < distgeo(1,j+1,2)
       if Dmax(1,j,1) < Dmax(1,j+1,2)
            Dmax(1,j,1) = Dmax(1,j+1,2);
       end
   end
   
   if distgeo(1,j,1) < distgeo(1,j-1,2)
       if Dmax(1,j,1) < Dmax(1,j-1,2)
            Dmax(1,j,1) = Dmax(1,j-1,2);
       end
   end
   
    if distgeo(1,j,1) < distgeo(2,j-1,2)
       if Dmax(1,j,1) < Dmax(2,j-1,2)
            Dmax(1,j,1) = Dmax(2,j-1,2);
       end
    end
   
     if distgeo(1,j,1) < distgeo(2,j+1,2)
       if Dmax(1,j,1) < Dmax(2,j+1,2)
            Dmax(1,j,1) = Dmax(2,j+1,2);
       end
     end
   
     if distgeo(1,j,1) < distgeo(2,j,2)
       if Dmax(1,j,1) < Dmax(2,j,2)
            Dmax(1,j,1) = Dmax(2,j,2);
       end
     end
   
 %%
 
     if distgeo(end,j,1) < distgeo(end,j+1,1)
       if Dmax(end,j,1) < Dmax(end,j+1,1)
            Dmax(end,j,1) = Dmax(end,j+1,1);
       end
   end
   
   if distgeo(end,j,1) < distgeo(end,j-1,1)
       if Dmax(end,j,1) < Dmax(end,j-1,1)
            Dmax(end,j,1) = Dmax(end,j-1,1);
       end
   end
   
    if distgeo(end,j,1) < distgeo(end-1,j-1,1)
       if Dmax(end,j,1) < Dmax(end-1,j-1,1)
            Dmax(end,j,1) = Dmax(end-1,j-1,1);
       end
    end
   
     if distgeo(end,j,1) < distgeo(end-1,j+1,1)
       if Dmax(end,j,1) < Dmax(end-1,j+1,1)
            Dmax(end,j,1) = Dmax(end-1,j+1,1);
       end
     end
   
     if distgeo(end,j,1) < distgeo(end-1,j,1)
       if Dmax(end,j,1) < Dmax(end-1,j,1)
            Dmax(end,j,1) = Dmax(end-1,j,1);
       end
     end
   
      if distgeo(end,j,1) < distgeo(end,j,2)
       if Dmax(end,j,1) < Dmax(end,j,2)
            Dmax(end,j,1) = Dmax(end,j,2);
       end
      end
   
      if distgeo(end,j,1) < distgeo(end,j+1,2)
       if Dmax(end,j,1) < Dmax(end,j+1,2)
            Dmax(end,j,1) = Dmax(end,j+1,2);
       end
   end
   
   if distgeo(end,j,1) < distgeo(end,j-1,2)
       if Dmax(end,j,1) < Dmax(end,j-1,2)
            Dmax(end,j,1) = Dmax(end,j-1,2);
       end
   end
   
    if distgeo(end,j,1) < distgeo(end-1,j-1,2)
       if Dmax(end,j,1) < Dmax(end-1,j-1,2)
            Dmax(end,j,1) = Dmax(end-1,j-1,2);
       end
    end
   
     if distgeo(end,j,1) < distgeo(end-1,j+1,2)
       if Dmax(end,j,1) < Dmax(end-1,j+1,2)
            Dmax(end,j,1) = Dmax(end-1,j+1,2);
       end
     end
   
     if distgeo(end,j,1) < distgeo(end-1,j,2)
       if Dmax(end,j,1) < Dmax(end-1,j,2)
            Dmax(end,j,1) = Dmax(end-1,j,2);
       end
     end
   
     %% change in z
     
        if distgeo(1,j,end) < distgeo(1,j+1,end)
       if Dmax(1,j,end) < Dmax(1,j+1,end)
            Dmax(1,j,end) = Dmax(1,j+1,end);
       end
   end
   
   if distgeo(1,j,end) < distgeo(1,j-1,end)
       if Dmax(1,j,end) < Dmax(1,j-1,end)
            Dmax(1,j,end) = Dmax(1,j-1,end);
       end
   end
   
    if distgeo(1,j,end) < distgeo(2,j-1,end)
       if Dmax(1,j,end) < Dmax(2,j-1,end)
            Dmax(1,j,end) = Dmax(2,j-1,end);
       end
    end
   
     if distgeo(1,j,end) < distgeo(2,j+1,end)
       if Dmax(1,j,end) < Dmax(2,j+1,end)
            Dmax(1,j,end) = Dmax(2,j+1,end);
       end
     end
   
     if distgeo(1,j,end) < distgeo(2,j,end)
       if Dmax(1,j,end) < Dmax(2,j,end)
            Dmax(1,j,end) = Dmax(2,j,end);
       end
     end
   
       if distgeo(1,j,end) < distgeo(1,j,end-1)
       if Dmax(1,j,end) < Dmax(1,j,end-1)
            Dmax(1,j,end) = Dmax(1,j,end-1);
       end
      end
      if distgeo(1,j,end) < distgeo(1,j+1,end-1)
       if Dmax(1,j,end) < Dmax(1,j+1,end-1)
            Dmax(1,j,end) = Dmax(1,j+1,end-1);
       end
   end
   
   if distgeo(1,j,end) < distgeo(1,j-1,end-1)
       if Dmax(1,j,end) < Dmax(1,j-1,end-1)
            Dmax(1,j,end) = Dmax(1,j-1,end-1);
       end
   end
   
    if distgeo(1,j,end) < distgeo(2,j-1,end-1)
       if Dmax(1,j,end) < Dmax(2,j-1,end-1)
            Dmax(1,j,end) = Dmax(2,j-1,end-1);
       end
    end
   
     if distgeo(1,j,end) < distgeo(2,j+1,end-1)
       if Dmax(1,j,end) < Dmax(2,j+1,end-1)
            Dmax(1,j,end) = Dmax(2,j+1,end-1);
       end
     end
   
     if distgeo(1,j,end) < distgeo(2,j,end-1)
       if Dmax(1,j,end) < Dmax(2,j,end-1)
            Dmax(1,j,end) = Dmax(2,j,end-1);
       end
     end
   
 %%
 
     if distgeo(end,j,end) < distgeo(end,j+1,end)
       if Dmax(end,j,end) < Dmax(end,j+1,end)
            Dmax(end,j,end) = Dmax(end,j+1,end);
       end
   end
   
   if distgeo(end,j,end) < distgeo(end,j-1,end)
       if Dmax(end,j,end) < Dmax(end,j-1,end)
            Dmax(end,j,end) = Dmax(end,j-1,end);
       end
   end
   
    if distgeo(end,j,end) < distgeo(end-1,j-1,end)
       if Dmax(end,j,end) < Dmax(end-1,j-1,end)
            Dmax(end,j,end) = Dmax(end-1,j-1,end);
       end
    end
   
     if distgeo(end,j,end) < distgeo(end-1,j+1,end)
       if Dmax(end,j,end) < Dmax(end-1,j+1,end)
            Dmax(end,j,end) = Dmax(end-1,j+1,end);
       end
     end
   
     if distgeo(end,j,end) < distgeo(end-1,j,end)
       if Dmax(end,j,end) < Dmax(end-1,j,end)
            Dmax(end,j,end) = Dmax(end-1,j,end);
       end
     end
   
      if distgeo(end,j,end) < distgeo(end,j,end-1)
       if Dmax(end,j,end) < Dmax(end,j,end-1)
            Dmax(end,j,end) = Dmax(end,j,end-1);
       end
      end
   
      if distgeo(end,j,end) < distgeo(end,j+1,end-1)
       if Dmax(end,j,end) < Dmax(end,j+1,end-1)
            Dmax(end,j,end) = Dmax(end,j+1,end-1);
       end
   end
   
   if distgeo(end,j,end) < distgeo(end,j-1,end-1)
       if Dmax(end,j,end) < Dmax(end,j-1,end-1)
            Dmax(end,j,end) = Dmax(end,j-1,end-1);
       end
   end
   
    if distgeo(end,j,end) < distgeo(end-1,j-1,end-1)
       if Dmax(end,j,end) < Dmax(end-1,j-1,end-1)
            Dmax(end,j,end) = Dmax(end-1,j-1,end-1);
       end
    end
   
     if distgeo(end,j,end) < distgeo(end-1,j+1,end-1)
       if Dmax(end,j,end) < Dmax(end-1,j+1,end-1)
            Dmax(end,j,end) = Dmax(end-1,j+1,end-1);
       end
     end
   
     if distgeo(end,j,end) < distgeo(end-1,j,end-1)
       if Dmax(end,j,end) < Dmax(end-1,j,end-1)
            Dmax(end,j,end) = Dmax(end-1,j,end-1);
       end
     end
  end
  end
     
  %check k edge

  for k=2:Az-1
  
     if distgeo(i,j,k) ~= 0
    
    if distgeo(1,1,k) < distgeo(1,1,k+1)
       if Dmax(1,1,k) < Dmax(1,1,k+1)
            Dmax(1,1,k) = Dmax(1,1,k+1);
       end
   end
   
   if distgeo(1,1,k) < distgeo(1,1,k-1)
       if Dmax(1,1,k) < Dmax(1,1,k-1)
            Dmax(1,1,k) = Dmax(1,1,k-1);
       end
   end
   
    if distgeo(1,1,k) < distgeo(1,2,k-1)
       if Dmax(1,1,k) < Dmax(1,2,k-1)
            Dmax(1,1,k) = Dmax(1,2,k-1);
       end
    end
   
     if distgeo(1,1,k) < distgeo(1,2,k+1)
       if Dmax(1,1,k) < Dmax(1,2,k+1)
            Dmax(1,1,k) = Dmax(1,2,k+1);
       end
     end
   
     if distgeo(1,1,k) < distgeo(1,2,k)
       if Dmax(1,1,k) < Dmax(1,2,k)
            Dmax(1,1,k) = Dmax(1,2,k);
       end
     end
   
       if distgeo(1,1,k) < distgeo(2,1,k)
       if Dmax(1,1,k) < Dmax(2,1,k)
            Dmax(1,1,k) = Dmax(2,1,k);
       end
      end
      if distgeo(1,1,k) < distgeo(2,1,k+1)
       if Dmax(1,1,k) < Dmax(2,1,k+1)
            Dmax(1,1,k) = Dmax(2,1,k+1);
       end
   end
   
   if distgeo(1,1,k) < distgeo(2,1,k-1)
       if Dmax(1,1,k) < Dmax(2,1,k-1)
            Dmax(1,1,k) = Dmax(2,1,k-1);
       end
   end
   
    if distgeo(1,1,k) < distgeo(2,2,k-1)
       if Dmax(1,1,k) < Dmax(2,2,k-1)
            Dmax(1,1,k) = Dmax(2,2,k-1);
       end
    end
   
     if distgeo(1,1,k) < distgeo(2,2,k+1)
       if Dmax(1,1,k) < Dmax(2,2,k+1)
            Dmax(1,1,k) = Dmax(2,2,k+1);
       end
     end
   
     if distgeo(1,1,k) < distgeo(2,2,k)
       if Dmax(1,1,k) < Dmax(2,2,k)
            Dmax(1,1,k) = Dmax(2,2,k);
       end
     end
   
 %%
 
     if distgeo(1,end,k) < distgeo(1,end,k+1)
       if Dmax(1,end,k) < Dmax(1,end,k+1)
            Dmax(1,end,k) = Dmax(1,end,k+1);
       end
   end
   
   if distgeo(1,end,k) < distgeo(1,end,k-1)
       if Dmax(1,end,k) < Dmax(1,end,k-1)
            Dmax(1,end,k) = Dmax(1,end,k-1);
       end
   end
   
    if distgeo(1,end,k) < distgeo(1,end-1,k-1)
       if Dmax(1,end,k) < Dmax(1,end-1,k-1)
            Dmax(1,end,k) = Dmax(1,end-1,k-1);
       end
    end
   
     if distgeo(1,end,k) < distgeo(1,end-1,k+1)
       if Dmax(1,end,k) < Dmax(1,end-1,k+1)
            Dmax(1,end,k) = Dmax(1,end-1,k+1);
       end
     end
   
     if distgeo(1,end,k) < distgeo(1,end-1,k)
       if Dmax(1,end,k) < Dmax(1,end-1,k)
            Dmax(1,end,k) = Dmax(1,end-1,k);
       end
     end
   
      if distgeo(1,end,k) < distgeo(2,end,k)
       if Dmax(1,end,k) < Dmax(2,end,k)
            Dmax(1,end,k) = Dmax(2,end,k);
       end
      end
   
      if distgeo(1,end,k) < distgeo(2,end,k+1)
       if Dmax(1,end,k) < Dmax(2,end,k+1)
            Dmax(1,end,k) = Dmax(2,end,k+1);
       end
   end
   
   if distgeo(1,end,k) < distgeo(2,end,k-1)
       if Dmax(1,end,k) < Dmax(2,end,k-1)
            Dmax(1,end,k) = Dmax(2,end,k-1);
       end
   end
   
    if distgeo(1,end,k) < distgeo(2,end-1,k-1)
       if Dmax(1,end,k) < Dmax(2,end-1,k-1)
            Dmax(1,end,k) = Dmax(2,end-1,k-1);
       end
    end
   
     if distgeo(1,end,k) < distgeo(2,end-1,k+1)
       if Dmax(1,end,k) < Dmax(2,end-1,k+1)
            Dmax(1,end,k) = Dmax(2,end-1,k+1);
       end
     end
   
     if distgeo(1,end,k) < distgeo(2,end-1,k)
       if Dmax(1,end,k) < Dmax(2,end-1,k)
            Dmax(1,end,k) = Dmax(2,end-1,k);
       end
     end
   
     %% change in z
     
      if distgeo(end,1,k) < distgeo(end,1,k+1)
       if Dmax(end,1,k) < Dmax(end,1,k+1)
            Dmax(end,1,k) = Dmax(end,1,k+1);
       end
   end
   
   if distgeo(end,1,k) < distgeo(end,1,k-1)
       if Dmax(end,1,k) < Dmax(end,1,k-1)
            Dmax(end,1,k) = Dmax(end,1,k-1);
       end
   end
   
    if distgeo(end,1,k) < distgeo(end,2,k-1)
       if Dmax(end,1,k) < Dmax(end,2,k-1)
            Dmax(end,1,k) = Dmax(end,2,k-1);
       end
    end
   
     if distgeo(end,1,k) < distgeo(end,2,k+1)
       if Dmax(end,1,k) < Dmax(end,2,k+1)
            Dmax(end,1,k) = Dmax(end,2,k+1);
       end
     end
   
     if distgeo(end,1,k) < distgeo(end,2,k)
       if Dmax(end,1,k) < Dmax(end,2,k)
            Dmax(end,1,k) = Dmax(end,2,k);
       end
     end
   
       if distgeo(end,1,k) < distgeo(end-1,1,k)
       if Dmax(end,1,k) < Dmax(end-1,1,k)
            Dmax(end,1,k) = Dmax(end-1,1,k);
       end
      end
      if distgeo(end,1,k) < distgeo(end-1,1,k+1)
       if Dmax(end,1,k) < Dmax(end-1,1,k+1)
            Dmax(end,1,k) = Dmax(end-1,1,k+1);
       end
   end
   
   if distgeo(end,1,k) < distgeo(end-1,1,k-1)
       if Dmax(end,1,k) < Dmax(end-1,1,k-1)
            Dmax(end,1,k) = Dmax(end-1,1,k-1);
       end
   end
   
    if distgeo(end,1,k) < distgeo(end-1,2,k-1)
       if Dmax(end,1,k) < Dmax(end-1,2,k-1)
            Dmax(end,1,k) = Dmax(end-1,2,k-1);
       end
    end
   
     if distgeo(end,1,k) < distgeo(end-1,2,k+1)
       if Dmax(end,1,k) < Dmax(end-1,2,k+1)
            Dmax(end,1,k) = Dmax(end-1,2,k+1);
       end
     end
   
     if distgeo(end,1,k) < distgeo(end-1,2,k)
       if Dmax(end,1,k) < Dmax(end-1,2,k)
            Dmax(end,1,k) = Dmax(end-1,2,k);
       end
     end
   
 %%
 
     if distgeo(end,end,k) < distgeo(end,end,k+1)
       if Dmax(end,end,k) < Dmax(end,end,k+1)
            Dmax(end,end,k) = Dmax(end,end,k+1);
       end
   end
   
   if distgeo(end,end,k) < distgeo(end,end,k-1)
       if Dmax(end,end,k) < Dmax(end,end,k-1)
            Dmax(end,end,k) = Dmax(end,end,k-1);
       end
   end
   
    if distgeo(end,end,k) < distgeo(end,end-1,k-1)
       if Dmax(end,end,k) < Dmax(end,end-1,k-1)
            Dmax(end,end,k) = Dmax(end,end-1,k-1);
       end
    end
   
     if distgeo(end,end,k) < distgeo(end,end-1,k+1)
       if Dmax(end,end,k) < Dmax(end,end-1,k+1)
            Dmax(end,end,k) = Dmax(end,end-1,k+1);
       end
     end
   
     if distgeo(end,end,k) < distgeo(end,end-1,k)
       if Dmax(end,end,k) < Dmax(end,end-1,k)
            Dmax(end,end,k) = Dmax(end,end-1,k);
       end
     end
   
      if distgeo(end,end,k) < distgeo(end-1,end,k)
       if Dmax(end,end,k) < Dmax(end-1,end,k)
            Dmax(end,end,k) = Dmax(end-1,end,k);
       end
      end
   
      if distgeo(end,end,k) < distgeo(end-1,end,k+1)
       if Dmax(end,end,k) < Dmax(end-1,end,k+1)
            Dmax(end,end,k) = Dmax(end-1,end,k+1);
       end
   end
   
   if distgeo(end,end,k) < distgeo(end-1,end,k-1)
       if Dmax(end,end,k) < Dmax(end-1,end,k-1)
            Dmax(end,end,k) = Dmax(end-1,end,k-1);
       end
   end
   
    if distgeo(end,end,k) < distgeo(end-1,end-1,k-1)
       if Dmax(end,end,k) < Dmax(end-1,end-1,k-1)
            Dmax(end,end,k) = Dmax(end-1,end-1,k-1);
       end
    end
   
     if distgeo(end,end,k) < distgeo(end-1,end-1,k+1)
       if Dmax(end,end,k) < Dmax(end-1,end-1,k+1)
            Dmax(end,end,k) = Dmax(end-1,end-1,k+1);
       end
     end
   
     if distgeo(end,end,k) < distgeo(end-1,end-1,k)
       if Dmax(end,end,k) < Dmax(end-1,end-1,k)
            Dmax(end,end,k) = Dmax(end-1,end-1,k);
       end
     end
   
     end
end
  




  end
     %%
  %checking Dmax
  
  
%   for i = 1:size(distgeo,1)-1;
%     
%     figure(100)
%     temp = Dmax(i,:,:);
%     imagesc(squeeze(temp));
%     axis equal
%     colorbar
%     drawnow
%     
%   end


%%
% calculate w
% rho = 1;
% mu = 1;
%   for i = 1:Ax
%       for j = 1:Ay
%          for k = 1:Az
%           
%              if distgeo(i,j,k) ~= 0
%                  
% 
%         w(i,j,k) =  Shape * R^2 * (rho/(8*mu)) * (2 * Dmax(i,j,k) * distgeo(i,j,k) - (distgeo(i,j,k))^2 );
% 
%              end
%           
%           end
%      end
%       end

     
%       clear distgeo
%       clear Dmax
      


disp('Dmax calculated')


end


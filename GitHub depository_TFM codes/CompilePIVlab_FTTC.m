
%FUNCTION: 
%1) to read and load into cell format - like what we get from PIVlab output.
%2) give an output.tif with overlayed velocity vectors on images, there are
%   Nframes*2 number of frames, because the 2(n-1)+1 frames will be the ref
%   frames, while the 2n frames will be each of the bead frames with cells.
%   


%%Note: 
%u(v)_TFM_xxx: are those without filtering of any forces
%however, the plotted in pictures are those filter based on the parameters
%below


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %Change parameters here
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       Nframes=32;
       Nstd_TFM=5;    %get rid of outliers from the mean for forces (based on force magnitude)  
       
       %take off n grid pts from border 
        ridNpts=0;
        jumppts=1;
        
        arrowhead_size=1;
        qscale=0.1;   %quiver scale, 0.01 for force
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %FINISH Change parameters here
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%


  x_TFM_GTest=cell(Nframes,1);
  y_TFM_GTest=cell(Nframes,1);
  u_TFM_GTest=cell(Nframes,1);   %force direction (x component)
  v_TFM_GTest=cell(Nframes,1);
   TF_MartielOutput_GTest=cell(Nframes,1);%%**
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%IMPORT Traction data TXT FILE SEQUENCES%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
      for i=1:Nframes
         fileName = strcat('Traction_PIV',num2str(i),'.txt');   %%**
         TF_MartielOutput_GTest{i,1}=importdata(fileName);
      end
      
      %%make it into PIVlab-style output
    
          xlist=sort(unique(TF_MartielOutput_GTest{1,1}(:,1)));
          ylist=sort(unique(TF_MartielOutput_GTest{1,1}(:,2)));
          dum=size(xlist);
          Numxpt=dum(1,1);
          dum=size(ylist);
          Numypt=dum(1,1);
          x_TFM_GTest{1,1}=[];
          y_TFM_GTest{1,1}=[];
          for j=1:Numypt
            x_TFM_GTest{1,1} =[x_TFM_GTest{1,1};xlist'];
          end
          for j=1:Numxpt
            y_TFM_GTest{1,1} =[y_TFM_GTest{1,1},ylist];
          end
          
          for i=2:Nframes
             x_TFM_GTest{i,1}=x_TFM_GTest{1,1}; 
             y_TFM_GTest{i,1}=y_TFM_GTest{1,1};
          end
          
     for i=1:Nframes
        for j=1:Numypt
            for k=1:Numxpt
                temp=TF_MartielOutput_GTest{i,1};
                temp=temp(temp(:,1)==x_TFM_GTest{i,1}(j,k),:);
                temp=temp(temp(:,2)==y_TFM_GTest{i,1}(j,k),:);
                
                u_TFM_GTest{i,1}(j,k)=temp(1,3);
                v_TFM_GTest{i,1}(j,k)=temp(1,4);
            end
        end
     end
    
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%FINISH import Traction data TXT FILE SEQUENCES%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%remove force outiers (how many std more than mean, using the force magnitude, for each frame)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
            %%this part gives NaN for outliers, and 0 for other values
                    Force_Magnitude=cell(Nframes,1);
                    for i=1:Nframes
                       Force_Magnitude{i,1}=sqrt(u_TFM_GTest{i,1}.*u_TFM_GTest{i,1}+v_TFM_GTest{i,1}.*v_TFM_GTest{i,1});
                       Force_Magnitude{i,1}(Force_Magnitude{i,1}>mean2(Force_Magnitude{i,1})+Nstd_TFM*std2(Force_Magnitude{i,1}))=NaN;
                       Force_Magnitude{i,1}=Force_Magnitude{i,1}*0;  %make values other than NaN to be 0;
                    end
            %%%%%%%%%%%%%%%%%%%%%%%%%

            %%%%set NaN for outliers of u, v based on force magnitude, and keep
            %%%%the others the same
                  u_TFM_filtered=cell(Nframes,1);
                  v_TFM_filtered=cell(Nframes,1);
                 for i=1:Nframes
                      u_TFM_filtered{i,1}= Force_Magnitude{i,1}+u_TFM_GTest{i,1};
                      v_TFM_filtered{i,1}= Force_Magnitude{i,1}+v_TFM_GTest{i,1}; 
                 end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%FINISH remove force outiers (how many std more than mean, using the force magnitude, for each frame)
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  
  
  
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 %%%%%%%%%%%%%%plot quiver vectors on images %%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if exist('output.tif', 'file') == 2
           delete 'output.tif' %this is to ensure there was no previous file with this name, so that we dont keep appending and exploding the original file end
end 

% Import map image
%srcFiles = dir('C:\Users\Thuan Beng\Desktop\test-nematicTFM\List\*.tif');  % the folder in which ur images exists
%length(srcFiles)
count=0;
for i = 1 : Nframes
       % filename = strcat('C:\Users\Thuan Beng\Desktop\test-nematicTFM\List\',srcFiles(i).name);
       
       for j=1:2
               %first part append the reference bead image
               count=count+1;
               if count<10
                    filename = strcat('C:\Users\uqjbroo4\Desktop\Output\ForPIV00',num2str(count),'.tif');   %%%%%%%%%%%%%%%%%%*********************************
               end
               if count>=10   && i <100
                    filename = strcat('C:\Users\uqjbroo4\Desktop\Output\ForPIV0',num2str(count),'.tif');   %%%%%%%%%%%%%%%%%%*********************************
               end
               if count>=100   && i <1000
                    filename = strcat('C:\Users\uqjbroo4\Desktop\Output\ForPIV',num2str(count),'.tif');   %%%%%%%%%%%%%%%%%%*********************************
               end
                I1 = imread(filename,1);
                f=figure('visible','off'); imshow(I1);
                hold on;

                x=x_TFM_GTest{i,1};
                y=y_TFM_GTest{i,1};
                u=u_TFM_filtered{i,1};
                v=v_TFM_filtered{i,1};



                dum=size(x);
                xrange=dum(1,2);
                yrange=dum(1,1);

                x_reduced=x(1+ridNpts:jumppts:yrange-ridNpts,1+ridNpts:jumppts:xrange-ridNpts);
                y_reduced=y(1+ridNpts:jumppts:yrange-ridNpts,1+ridNpts:jumppts:xrange-ridNpts);
                u_reduced=u(1+ridNpts:jumppts:yrange-ridNpts,1+ridNpts:jumppts:xrange-ridNpts);
                v_reduced=v(1+ridNpts:jumppts:yrange-ridNpts,1+ridNpts:jumppts:xrange-ridNpts);


                h=quiver(x_reduced,y_reduced,u_reduced,v_reduced,0,'c');

                hU = get(h,'UData') ;
                hV = get(h,'VData') ;
                set(h,'UData',qscale*hU,'VData',qscale*hV)

                %adjust_quiver_arrowhead_size(h, arrowhead_size);   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                f = getframe(gca);
                im = frame2im(f);


                outputFileName = 'output.tif';
                imwrite(im, outputFileName, 'WriteMode', 'append');
                  close
       end


%         %second part append the cells image, on the ith frame
% 
%         I2 = imread(filename,2);
%         f=figure('visible','off'); imshow(I2);
%         hold on;
% 
%         h=quiver(x_reduced,y_reduced,u_reduced,v_reduced,0,'c');
%         qscale=0.01;
%         hU = get(h,'UData') ;
%         hV = get(h,'VData') ;
%         set(h,'UData',qscale*hU,'VData',qscale*hV)
% 
%         adjust_quiver_arrowhead_size(h, arrowhead_size);
%         f = getframe(gca);
% 
%         im = frame2im(f);
%         imwrite(im, outputFileName, 'WriteMode', 'append');
%           close
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 %%%%%%%%%%%%%%FINISH plot quiver vectors on images %%%%%%%%%%%%%%%%%%%%%%%%%%%% 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  



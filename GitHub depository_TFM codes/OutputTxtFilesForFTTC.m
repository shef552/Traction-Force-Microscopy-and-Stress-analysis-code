%convert PIVlab output into an output txtfiles recognized by FTTC plugin in Image J
%to get TFM. FTTC reads displacment outputs in microns.

%Input: (PIVlab output
%x: cell format. 
%y:
%u_filtered
%v_filtered.
%Conv_pixtoum: conversion factor to change pixels to microns  (for 20X imq
%movie, it is 0.32 um/pix).

%No output in matlab variable format, but txtfiles are generated. 

%Note: Becareful of output units!!!

function OutputTxtFilesForFTTC(x,y,u_filtered,v_filtered,Conv_pixtoum)
    
   dum=size(x);
   Numframes=dum(1,1);
   dum=size(x{1,1});
   xrange=dum(1,2);
   yrange=dum(1,1);
   
   for i=1:Numframes
       temp=[];
       for j=1:yrange
           for k=1:xrange
               temp=[temp;[x{i,1}(j,k),y{i,1}(j,k),u_filtered{i,1}(j,k),v_filtered{i,1}(j,k)]];
           end
       end
       temp(:,3:4)=temp(:,3:4)*Conv_pixtoum;
       
       fileID = fopen(['PIV',num2str(i),'.txt'],'w');
       fprintf(fileID,'%6.0f %6.0f %6.12f %6.12f\n',temp');  %important to be temp'
       fclose(fileID);
   end
   
   
end
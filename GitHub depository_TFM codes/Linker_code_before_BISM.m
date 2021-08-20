for i=1:Nframes
    f=sprintf('frame%d', i);
    traction.(f)=struct('tx',u_TFM_GTest{i},'ty',v_TFM_GTest{i});
    
end

save('Traction_field.mat','traction')


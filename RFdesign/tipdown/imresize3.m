function  output = imresize3(input,scale) 
for iz = 1:size(input,3)
   output(:,:,iz) = imresize(input(:,:,iz),scale);
end 
   
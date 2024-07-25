function slopes = fieldToSlopes(field,lenslets_across,lensletMask,targetTotalPower,rn)

%Shack-Hartment wavefront sensor function to measured slopes from the shift
%in intensity peaks of each lenslet 

% input field is a 192x192 matrix. If each lenslet is 16x16 accross, the
% field only covers 12 lenslets across. 

lensletSizePx = size(field)/lenslets_across; %number of lenslets that covers the field.

if ~mod(lensletSizePx,2) %If # lenslets across is even, pad paddedLenslet array by 1 for fft2 
    pad = lensletSizePx(1) +1; 
else 
    pad = lensletSizePx(1); 
end 

if ~exist('targetTotalPower','var')
    targetTotalPower = 0;
end

slopes = zeros(sum(sum(lensletMask)),2);
validLensletCounter = 0;
% disp(targetTotalPower);
%return 
for lX = 1:lenslets_across
    for lY = 1:lenslets_across
%         if lY == 8
%             'boo'
%         end
        currLenslet = field((1:lensletSizePx(1)) + lensletSizePx(1)*(lY-1),(1:lensletSizePx(2)) + lensletSizePx(2)*(lX-1));
        
        paddedLenslet = currLenslet;
        

        if ~lensletMask(lX,lY)
            continue;
        end

        validLensletCounter = validLensletCounter +1;
        focusedField = fftshift(fft2(paddedLenslet,pad,pad));
        
        intensity = abs(focusedField).^2; %very important!!!!
        totalPower = sum(sum(intensity));
        % noise application block
        if targetTotalPower > 0
            % add poiss shot noise
            %intensity = poissrnd(intensity*targetTotalPower/totalPower);
            % add read noise with zero mean and std rn
            intensity = intensity + rn.*randn(size(focusedField));
            % clip to 0
            intensity = max(intensity,0);
            % clip to only values above 2 sigma of read noise (is this ideal?)
            intensity = intensity - 2*rn;
            intensity = max(intensity,0);
        end
        
        % centroiding
        [X,Y] = meshgrid(1:size(intensity,1),1:size(intensity,2));
        X = X - floor(size(intensity,2)/2)-1;
        Y = Y - floor(size(intensity,1)/2)-1;
        
        cenX = sum(sum(X.*intensity))/sum(intensity(:));
        cenY = sum(sum(Y.*intensity))/sum(intensity(:));
         %figure(1); imagesc(intensity); colormap gray; [cenX cenY]
         %pause(0.1)
         %slopes((0:1)+validLensletCounter) = [cenX,cenY];
        slopes(validLensletCounter,1) = cenX;
        slopes(validLensletCounter,2) = cenY;
        
        %csvwrite(sprintf('./receiverdata/intensity/test_X%i_Y%i.csv',lX, lY),intensity);
        
%         FID = fopen(sprintf('./receiverdata/intensity/test_int_X%i_Y%i.bin',lX, lY), 'w+'); 
%         fwrite(FID, intensity, 'float32', 'b'); fclose(FID);
%         
%         FID = fopen(sprintf('./receiverdata/center/test_center_X%i_Y%i.bin',lX, lY), 'w+'); 
%         fwrite(FID, [cenX,cenY], 'float32', 'b'); fclose(FID);
    end
end



% slopes = slopes(1:(validLensletCounter*2));
function [minVal,minId] = getNearbySigma(sizeWISEminSigma,k)
% Get the better length for a signal
% It has been changed to EMD implementation
% Now we are getting into a very minmum value since we are working on distance (by EMD) instead of similarity value 
% Dated: 8th August 2019
% by munna
MinSigmas = round(sizeWISEminSigma,4);     % round to 4 decimal points 
[uMinSigmas,uIndx] = unique(MinSigmas);    % Take unique entries where we always get lower index for a identical simScore 
[vals,mIndx]       = mink(uMinSigmas,k);   % 5 MinSigmas
idss = uIndx(mIndx);                       % Real index 
%[vals,idss]= maxk(sizeWISEmaxSigma,k); 
minId  = idss(1);
minVal = vals(1);

for i =1:length(idss)
    if(idss(i)<minId)
       
        % Take lower index if the simScore is equal
        %if(vals(i)<=minVal)
         %   minId = idss(i);
        %end
        %else
        % ------------------------
        gainLen = ((minId-idss(i))*100)/minId; 
        increasedDist = ((vals(i)-minVal)*100)/minVal;
        
        if((vals(i)-vals(1))<= 5 && increasedDist<=2 && gainLen>=50)
            % For 20 % loss in sim or increase in distance
            % Gain in index > 50% which means almost 14 index earlier
            % We have total 29 lengths 22-50
            % So, 14 is almost half of the curve
            minId = idss(i);
            minVal = vals(i);
        end
    end
end
end
function [miNoLight_allCells,miMidLight_allCells,miHighLight_allCells] = load_ResponseData(responses_expt)
    response_noLightOrig = responses_expt([1,4,7,10,13,16,19],:);
    response_midLightOrig = responses_expt([2,5,8,11,14,17,20],:);
    response_highLightOrig = responses_expt([3,6,9,12,15,18,21],:);
    response_noLightNorm = zeros(size(response_noLightOrig));
    response_midLightNorm = zeros(size(response_midLightOrig));
    response_highLightNorm = zeros(size(response_highLightOrig));
    
    numSoundResponsiveCells = size(responses_expt,2);
    miNoLight_allCells = zeros(numSoundResponsiveCells,1);
    miMidLight_allCells = zeros(numSoundResponsiveCells,1);
    miHighLight_allCells = zeros(numSoundResponsiveCells,1);
    for iCell = 1:numSoundResponsiveCells
        response_noLightNorm(:,iCell) = response_noLightOrig(:,iCell)/max(abs(response_noLightOrig(:,iCell)));
        response_midLightNorm(:,iCell) = response_midLightOrig(:,iCell)/max(abs(response_midLightOrig(:,iCell)));
        response_highLightNorm(:,iCell) = response_highLightOrig(:,iCell)/max(abs(response_highLightOrig(:,iCell)));
        
        miNoLight_allCells(iCell) = monotonicity_rawData(response_noLightNorm(:,iCell));
        miMidLight_allCells(iCell) = monotonicity_rawData(response_midLightNorm(:,iCell));
        miHighLight_allCells(iCell) = monotonicity_rawData(response_highLightNorm(:,iCell));
    end
end
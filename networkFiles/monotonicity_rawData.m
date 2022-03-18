function mi = monotonicity_rawData(response)
    %Using the monotonicity index of paper - Level-Tuned Neurons in Primary Auditory Cortex Adapt Differently to Loud versus Soft Sounds
    mi = abs(response(end)-response(1))/abs(max(response)-min(response));
end
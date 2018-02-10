function compileInference(package)

  if nargin==0, 
    package={'inference', 'full_factor_gbp', 'factor_gbp', 'gbp_preprocess', 'gbp'};
  end
  
  % compile c_inference - the whole inference package
  if any(strcmp('inference', package)),
    mex -O c_inference.cpp fillMethods.cpp InferenceAlgorithm.cpp Loopy.cpp LogLoopy.cpp LoopySTime.cpp LogLoopySTime.cpp PairsGBP.cpp LogPairsGBP.cpp MeanField.cpp LogMeanField.cpp GBP.cpp LogGBP.cpp MonteCarlo.cpp Gibbs.cpp Wolff.cpp SwendsenWang.cpp Metropolis.cpp GBPPreProcessor.cpp Region.cpp RegionLevel.cpp MRF.cpp PottsMRF.cpp
  end
  
  % compile c_factor_gbp
  if any(strcmp('factor_gbp', package)),
    mex -O c_factor_gbp.cpp fillMethods.cpp InferenceAlgorithm.cpp GBP.cpp Region.cpp RegionLevel.cpp MRF.cpp PottsMRF.cpp
  end
  % compile c_full_factor_gbp
  if any(strcmp('full_factor_gbp', package)),
    mex -O c_full_factor_gbp.cpp fillMethods.cpp InferenceAlgorithm.cpp GBP.cpp GBPPreProcessor.cpp Region.cpp RegionLevel.cpp MRF.cpp PottsMRF.cpp
  end
  % compile c_gbp_preprocess - makes only the preprocessing for gbp
  if any(strcmp('gbp_preprocess', package)),
    mex -O c_gbp_preprocess.cpp fillMethods.cpp GBPPreProcessor.cpp Region.cpp RegionLevel.cpp MRF.cpp PottsMRF.cpp
  end

  % compile c_gbp - makes only the inference for gbp, using data resulted by c_gbp_preprocess
  if any(strcmp('gbp', package)),
    mex -O c_gbp.cpp fillMethods.cpp InferenceAlgorithm.cpp GBP.cpp Region.cpp RegionLevel.cpp MRF.cpp PottsMRF.cpp
  end

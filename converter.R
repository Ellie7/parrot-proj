install.packages("matconv")
library(matconv)

matIn <- c(%CorrRates: this is a program to generate sets of correlated
            % vital rates.
            clear all;
            %*****************Simulation Parameters***********************
              % the basic information on the vital rates first is for a beta, second
            %for a lognormal, and then a strectched beta.
            vrmeans= [0.77,2,5]; %means
            vrvars= [0.02,0.5,2]; %variances -- not standard deviations
            % minimum and maximum values for each vital rate: zeros are place holders
            %for rates that are not stretched betas ;
            vrmins= [0 0 0];
            vrmaxs=[0 0 24];
            %then a useable correlation matrix (corrected if original was not good)
            corrmx=[...
                    1 0.3 -0.2;
                    0.3 1 -0.6;
                    -0.2 -0.6 1]
            tmax = 200; %number of years of vital rates to simulate;
            %*************************************************************
              np = length(vrmeans);
            results = [];
            randn('state',sum(100*clock)); %this seeds the random number generator
            %find the eigenvalues, eee, and eigenvectors of the correlations
            [uuu,eee] = eig(corrmx);
            %Calculate z12, the matrix to use to make uncorrelated st. normals correlated
            z12 = uuu*(sqrt(abs(eee)))*uuu;
            for tt=1:tmax
            normvals = randn(np,1); %make a set of random st. normal values.
            corrnorms = z12*normvals; %make correlated st. normals
            %get the beta vital rate with the same Fx as the first normal:
            vrate1 = betaval(vrmeans(1),(vrvars(1))^0.5,stnormfx(corrnorms(1)));
            %convert the second norm into a corresponding log normal:
            vrate2 = lnorms(vrmeans(2),vrvars(2),corrnorms(2));
            %get the stretched beta vital rate with the same Fx as the 3rd normal:
            vrate3 = stretchbetaval(vrmeans(3),(vrvars(3))^0.5,...
            vrmins(3), vrmaxs(3),stnormfx(corrnorms(3)));
            results = [results;vrate1,vrate2,vrate3];
            end; %tt
            meanrates = mean(results)
            variances = var(results)
            correlations = corrcoef(results))

mat2r(matIn, verbose = 0)$rCode

# "xlsReadPretty <- function(...){" 
# "\tvarargin <- list(...)"
# "  didThing <- 1*3"
# "  dat <- didThing / 3"
#"\treturn(dat)"
#"}"
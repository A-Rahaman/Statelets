function [Sdyn, dSdyn] = create_surrogate_FNCdyn(Fdyn, METHOD)
% FNCdyn should be [nSub X nWin x nFeatures]

[M, Nwin, nCC] = size(Fdyn);

Sdyn = zeros(size(Fdyn));

for sub = 1:M
    fprintf('Working on subject %d of %d\n', sub, M)
    s = squeeze(Fdyn(sub,:,:));
    if strcmp(METHOD, 'consistent')
        ms = repmat(mean(s), [Nwin 1]);
        sperm = Surrogate_FourierTransform(s-ms);
        % randomize the phase consistently in the fourier domain
        sperm = sperm + ms;
        Sdyn(sub,:,:) = sperm;
    elseif strcmp(METHOD, 'inconsistent')
        for ii = 1:nCC
            si = s(:,ii);
            ms = mean(si);
            sperm = Surrogate_FourierTransform(si-ms);
            sperm = sperm + ms;
            Sdyn(sub,:,ii) = sperm;
        end
    else
        disp('Method not recognized')
    end           
end

% remove the mean for clustering
dSdyn = Sdyn - repmat(mean(Sdyn,2), [1, Nwin, 1]);
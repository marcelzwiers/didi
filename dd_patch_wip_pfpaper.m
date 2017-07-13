%% PF artefact paper

% Load the phase-data (only once as it is persistent in dd_patch)
if isempty(PhaseVol)
    FNames   = handles.input(handles.DWISel,:);
    PathStr  = fileparts(FNames(1,:));
    FNames   = dir([PathStr filesep 'phase' filesep '*.img']);
    FNames   = [repmat([PathStr filesep 'phase' filesep], size(FNames)) strvcat(FNames.name)];
    PhaseVol = permute(spm_read_vols(spm_vol(FNames(handles.DWISel,:))), handles.tzyx);
    PhaseVol = pi*PhaseVol/4096;  % [-4096, 4096] => [-pi, pi]
end

% Find the corresponding repetitions
Qi  = find(sum(handles.q==repmat(BVec, [size(handles.q,1) 1]),2));
nQi = numel(Qi);
% Qi([1 find(Qi==qi)]) = Qi([find(Qi==qi) 1]);   % Put the current slice at the front

% Reconstruct the k-space images ([y,x,n])
for n = 1:nQi
    DWI(:,:,n)         = shiftdim(handles.Vol(Qi(n), zi, :, :));
    PhaseIm(:,:,n)     = shiftdim(PhaseVol(Qi(n), zi, :, :));
    %[dxPhase dyPhase]  = gradient(PhaseIm(:,:,n)); % NB: No unwrapping
    %PhaseGradIm(:,:,n) = sign(dyPhase) .* sqrt(dxPhase.^2 + dyPhase.^2); % with sign
    [dum PhaseGradIm(:,:,n)] = gradient(unwrap(PhaseIm(:,:,n))); % Unwrapped and derivated along the PF/PE-direction
    KSpace(:,:,n)      = fft2(DWI(:,:,n) .* exp(1i*PhaseIm(:,:,n)));
    KSpaceMag          = ifftshift(imfilter(fftshift(abs(KSpace(:,:,n))), ...
                         fspecial('gaussian',6,1)));
    KSpaceProf(:,:,n)  = fftshift([KSpaceMag(:,1) KSpaceMag(1,:)'], 1); % Powerprofile along the axes
    KSpace(:,:,n)      = fftshift(KSpace(:,:,n));
end

% Clean up the images
PhaseImUnwrap = PhaseIm;
for n = 1:nQi
    [PhaseImUnwrap(:,:,n) LowResMagIm(:,:,n) LowResPhaseIm(:,:,n)] = ...
        unwrap_homodyne(DWI(:,:,n), PhaseIm(:,:,n), floor(min(size(DWI(:,:,n)))/8));  % Remove the low-frequency wraps
end
%PhaseGradIm = smooth3(PhaseGradIm, 'gaussian', [3 3 1]);

%% Plot the phase and k-space data
figure(1)
delete(findall(gcf, 'Tag', 'SliceSelArrow'));
colormap('gray')
set(gcf, 'Position', [5 52 954 797])
NCol = 4;

subplot(1, NCol, 1)
    DWI = permute(DWI(YRange, XRange, :), [1 3 2]); % Plot all directions in 1 montage
    imagesc(reshape(DWI, [size(DWI,1)*nQi size(DWI,3)]))
    Pos = get(gca, 'Position');
    X = [0.4*Pos(1) 0.6*Pos(1)];
    Y = (Pos(2) + Pos(4)*(find(flipud(Qi)==qi)-0.5)/nQi) * [1 1];
    annotation(gcf, 'arrow', X, Y, 'Tag', 'SliceSelArrow');
    %colorbar
    axis image off
    title('Magnitude')
    hold on
    for n = 1:nQi     % Plot the artefact edge
        for m = 1:numel(EdgeX)
            plot(EdgeY{m}, EdgeX{m}+(n-1)*size(DWI,1), '-c')
        end
    end
    hold off

subplot(1, NCol, 2)
    PhaseIm = permute(PhaseIm(YRange, XRange, :), [1 3 2]);
    imagesc(reshape(PhaseIm, [size(PhaseIm,1)*nQi size(PhaseIm,3)]))
    %colorbar
    axis image off
    title('Phase')
    hold on
    for n = 1:nQi     % Plot the artefact edge
        for m = 1:numel(EdgeX)
            plot(EdgeY{m}, EdgeX{m}+(n-1)*size(DWI,1), '-c')
        end
    end
    hold off
    
subplot(1, NCol, 3)
    PhaseGradIm = permute(PhaseGradIm(YRange, XRange, :), [1 3 2]);
    imagesc(reshape(PhaseGradIm, [size(PhaseGradIm,1)*nQi size(PhaseGradIm,3)]))
    %colorbar
    axis image off
    title('\partialPhase/\partialy')
    hold on
    for n = 1:nQi     % Plot the artefact edge
        for m = 1:numel(EdgeX)
            plot(EdgeY{m}, EdgeX{m}+(n-1)*size(DWI,1), '-c')
        end
    end
    hold off

subplot(1, NCol, 4)
    KSpaceFlip = permute(KSpace, [1 3 2]);            % => [y',n,x]
    KSpaceFlip = reshape(KSpaceFlip, [size(KSpaceFlip,1)*nQi size(KSpaceFlip,3)]); % => [y'n,x]
    KSpace     = permute(flipdim(KSpace,1), [1 3 2]); % => [y,n,x]
    KSpace     = reshape(KSpace, [size(KSpace,1)*nQi size(KSpace,3)]); % => [yn,x]
    imagesc(log(abs(KSpace)))
    %colorbar
    axis image off
    title('log(|k|)')

% subplot(1, NCol, 5)
%     imagesc(log(abs(KSpace-conj(KSpaceFlip))))
%     %colorbar
%     axis image off
%     title('log(|k-conj(-k)|)')

    
%% Plot the low-frequency information
figure(2)
delete(findall(gcf, 'Tag', 'SliceSelArrow'));
colormap('gray')
%set(gcf, 'Position', [5 52 954 797])
NCol = 3;

%subplot(1, NCol, 1)
%     PhaseImUnwrap = permute(PhaseImUnwrap(YRange, XRange, :), [1 3 2]);
%     imagesc(reshape(PhaseImUnwrap, [size(PhaseImUnwrap,1)*nQi size(PhaseImUnwrap,3)]))
%     %colorbar
%     axis image off
%     title('Unwrapped phase')
%     hold on
%     for n = 1:nQi     % Plot the artefact edge
%         for m = 1:numel(EdgeX)
%             plot(EdgeY{m}, EdgeX{m}+(n-1)*size(DWI,1), '-c')
%         end
%     end
%     hold off

subplot(1, NCol, 2)
    LowResMagIm = permute(LowResMagIm(YRange, XRange, :), [1 3 2]);
    imagesc(reshape(LowResMagIm, [size(LowResMagIm,1)*nQi size(LowResMagIm,3)]))
    %colorbar
    axis image off
    title('Magnitude')
    hold on
    for n = 1:nQi     % Plot the artefact edge
        for m = 1:numel(EdgeX)
            plot(EdgeY{m}, EdgeX{m}+(n-1)*size(DWI,1), '-c')
        end
    end
    hold off

subplot(1, NCol, 3)
    LowResPhaseIm = permute(LowResPhaseIm(YRange, XRange, :), [1 3 2]);
    imagesc(reshape(LowResPhaseIm, [size(LowResPhaseIm,1)*nQi size(LowResPhaseIm,3)]))
    %colorbar
    axis image off
    title('Phase')
    hold on
    for n = 1:nQi     % Plot the artefact edge
        for m = 1:numel(EdgeX)
            plot(EdgeY{m}, EdgeX{m}+(n-1)*size(DWI,1), '-c')
        end
    end
    hold off

    
%% Plot the FWHM of k-space
figure(3), clf
NCol = 4;

KSpaceProfAv = squeeze(mean([KSpaceProf flipdim(KSpaceProf(:,2,:),1)], 2));

subplot(NCol, 1, 1)
    plot(KSpaceProf(:,:,1))
    hold on
    plot(KSpaceProfAv(:,1), 'LineWidth',2)
    set(gca, 'YScale', 'Log')
    vertlijn(ceil(size(KSpaceProf,1)/2), '--')
    vertlijn(ceil(size(KSpaceProf,1)*[5/8 3/8]), ':')
    ylabel('|k|')
    legend('y-dir', 'x-dir', 'mean')

subplot(NCol, 1, 2)
    plot(KSpaceProf(:,:,2))
    hold on
    plot(KSpaceProfAv(:,2), 'LineWidth',2)
    set(gca, 'YScale', 'Log')
    vertlijn(ceil(size(KSpaceProf,1)/2), '--')
    vertlijn(ceil(size(KSpaceProf,1)*[5/8 3/8]), ':')
    ylabel('|k|')

subplot(NCol, 1, 3)
    plot(KSpaceProf(:,:,3))
    hold on
    plot(KSpaceProfAv(:,3), 'LineWidth',2)
    set(gca, 'YScale', 'Log')
    vertlijn(ceil(size(KSpaceProf,1)/2), '--')
    vertlijn(ceil(size(KSpaceProf,1)*[5/8 3/8]), ':')
    ylabel('|k|')
    
subplot(NCol, 1, 4)
    plot(KSpaceProfAv)
    set(gca, 'YScale', 'Log')
    vertlijn(ceil(size(KSpaceProf,1)/2), '--')
    vertlijn(ceil(size(KSpaceProf,1)*[5/8 3/8]), ':')
    ylabel('|k|')
    legend('M1', 'M2', 'M3')

% figure(4), clf
%     plot(mean(KSpaceProfAvy))
%     set(gca, 'YScale', 'Log')
%     vertlijn(ceil(size(KSpaceProf,1)/2), '--')
%     vertlijn(ceil(size(KSpaceProf,1)*[5/8 3/8]), ':')
%     ylabel('|k|')
%     legend('M1', 'M2', 'M3')


return
figure(5),clf % ISMRM abstract

%% Plot the phase and k-space data
colormap('gray')
set(gcf, 'Position', [5 52 954 797])

Br = 0.13;
axes('Position', [0.1 0.1 Br 0.8])
    imagesc(reshape(DWI, [size(DWI,1)*nQi size(DWI,3)]))
    axis image off
    title('Magnitude')
    hold on
    for n = 1:nQi     % Plot the artefact edge
        for m = 1:numel(EdgeX)
            plot(EdgeY{m}, EdgeX{m}+(n-1)*size(DWI,1), '-c')
        end
    end
    hold off

axes('Position', [0.1+Br 0.1 Br 0.8])
    imagesc(reshape(PhaseIm, [size(PhaseIm,1)*nQi size(PhaseIm,3)]))
    %colorbar
    axis image off
    title('Phase')
    hold on
    for n = 1:nQi     % Plot the artefact edge
        for m = 1:numel(EdgeX)
            plot(EdgeY{m}, EdgeX{m}+(n-1)*size(DWI,1), '-c')
        end
    end
    hold off
    
axes('Position', [0.1+2*Br 0.1 Br 0.8])
    imagesc(reshape(PhaseGradIm, [size(PhaseGradIm,1)*nQi size(PhaseGradIm,3)]))
    %colorbar
    axis image off
    title('\partialPhase/\partialy')
    hold on
    for n = 1:nQi     % Plot the artefact edge
        for m = 1:numel(EdgeX)
            plot(EdgeY{m}, EdgeX{m}+(n-1)*size(DWI,1), '-c')
        end
    end
    hold off

axes('Position', [0.1+3*Br 0.1 1.226*Br 0.8])
    imagesc(log(abs(KSpace)))
    %colorbar
    axis image off
    title('log(|k|)')

    
%% Plot the low-frequency information

axes('Position', [0.1+4.2*Br 0.1 Br 0.8])
    imagesc(reshape(LowResMagIm, [size(LowResMagIm,1)*nQi size(LowResMagIm,3)]))
    %colorbar
    axis image off
    title('Magnitude')
    hold on
    for n = 1:nQi     % Plot the artefact edge
        for m = 1:numel(EdgeX)
            plot(EdgeY{m}, EdgeX{m}+(n-1)*size(DWI,1), '-c')
        end
    end
    hold off

axes('Position', [0.1+5.2*Br 0.1 Br 0.8])
    imagesc(reshape(LowResPhaseIm, [size(LowResPhaseIm,1)*nQi size(LowResPhaseIm,3)]))
    %colorbar
    axis image off
    title('Phase')
    hold on
    for n = 1:nQi     % Plot the artefact edge
        for m = 1:numel(EdgeX)
            plot(EdgeY{m}, EdgeX{m}+(n-1)*size(DWI,1), '-c')
        end
    end
    hold off

function [roi_points, intens] = clicky(movie_in, refimg);
% function [roi_points, intens] = clicky(movie_in, refimg);
% movie_in: input movie.
% refimg: optional image (black and white or color) on which to click to
% select ROIs.
% Click to define an ROI, right click to complete.  Plots the intensity
% trace averaged over the ROI.  Keeps generating ROIs until you generate an
% ROI with no points (i.e. just a right click on the image).
% Returns coordinates of all ROIs and the intensity traces.
% AEC and JMK, 12 Oct. 2011.

if isa(movie_in,'vm')
    movie_in = movie_in.toimg.data;
end

if nargin < 2;
    refimg = mean(movie_in, 3);
end;

movie_in = tovec(movie_in);

figure
subplot(1,3,1)
imshow(refimg, [], 'InitialMagnification', 'fit')
hold on;

[ysize, xsize] = size(refimg(:,:,1)); 
npts = 1;
colorindex = 0;
order = get(gca,'ColorOrder');
nroi = 1;
intens = [];
  [x, y] = meshgrid(1:xsize, 1:ysize);
 while(npts > 0)

  subplot(1,3,1)
  [xv, yv] = (getline(gca, 'closed'));
  if size(xv,1) < 3  % exit loop if only a line is drawn
  break
  end
  inpoly = inpolygon(x,y,xv,yv);
  
 %draw the bounding polygons and label them
  currcolor = order(1+mod(colorindex,size(order,1)),:);
  plot(xv, yv, 'Linewidth', 1,'Color',currcolor);
  text(mean(xv),mean(yv),num2str(colorindex+1),'Color',currcolor,'FontSize',12);
  
%   itrace = squeeze(sum(sum(movie_in.*repmat(inpoly, [1, 1, nframes]))))/sum(inpoly(:));
%   itrace = tovec(movie_in)'*inpoly(:)/nnz(inpoly); % 2016 VP
%   itrace = mean(movie_in(inpoly,:))'; % 2017 VP has bug for 1-pixel ROIs
    itrace = mean(movie_in(inpoly,:),1)'; % 2019 VP

  subplot(1,3,2:3) % plot the trace
  hold on;
  plot(itrace,'Color',currcolor);
  colorindex = colorindex+1;
  
  intens = [intens; itrace'];
  roi_points{nroi} = [xv, yv];
  nroi = nroi + 1;
 end
 intens = intens';
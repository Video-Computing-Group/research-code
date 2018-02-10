function [cordnew1]=getcentroid_f(phi,alpha,cellarea_threshold)

[ny,nx]=size(phi);

for i=1:ny
    for j=1:nx
        if phi(i,j)==1
            JJ(i,j)=255;        
        end
    end
end
%colormap([1 0 0])
%imshow(J);
%J= bwmorph(JJ, 'thin' , Inf);
%figure, imshow(J);
% BW = edge(JJ,'canny',[],alpha,'nothinning'); %Edge Transformation

%figure,imshow(BW);


cord=regionprops(bwlabeln(phi),'pixellist','PixelIdxList', 'BoundingBox','area','ConvexArea','EulerNumber','ConvexHull','Extrema','centroid','Solidity','Extent','EquivDiamete','Perimeter','MajorAxisLength', 'MinorAxisLength','Orientation','FilledArea');

m=length(cord);
j=1;


for ii=1:m
    for jj=1:m
        graph_matrix(ii,jj)=sqrt((cord(ii).Centroid(1)-cord(jj).Centroid(1))^2+(cord(ii).Centroid(2)-cord(jj).Centroid(2))^2);
    end
end

for iii=1:m
   temp=sort(graph_matrix(iii,:));
   lengthsum(iii)=temp(1)+temp(2)+temp(3)+temp(4)+temp(5)+temp(6);
end
   
    
for i=1:m
    if (cord(i).BoundingBox(3)>cellarea_threshold)&&(cord(i).BoundingBox(4)>cellarea_threshold)&&(cord(i).BoundingBox(3)<100)&&(cord(i).BoundingBox(4)<100)&&(lengthsum(i)<300)&&(cord(i).MinorAxisLength>10)
        cordnew(j)=cord(i);
        j=j+1;
    else
        j=j;
    end
end

m_new=length(cordnew);
for i=1:m_new
    x(i)=cordnew(i).Centroid(1);
    y(i)=cordnew(i).Centroid(2);
end

peripheral=convhull(x,y,{'Qt','Pp'});

j=1;
for i=1:m_new
    aa=find(peripheral==i);
    if isempty(aa)==1
        cordnew2(j)=cordnew(i);
        j=j+1;
    else
        j=j;
    end
end

imshow(phi); hold on;             %% Changed by Anirban 2/21/11

n=length(cordnew);

cordnew1=struct('Area', {0}, 'Centroid', {0}, 'BoundingBox', {0}, 'MajorAxisLength', {0}, 'MinorAxisLength', {0},'Orientation', {0}, ...
    'ConvexHull', {0}, 'ConvexArea', {0}, 'FilledArea',  {0}, 'EulerNumber', {0}, 'Extrema', {0}, 'EquivDiameter', {0}, 'Solidity', {0},...
    'Extent', {0}, 'PixelIdxList', {0}, 'PixelList', {0}, 'Perimeter', {0}, 'SliceNo', {0}, 'ParentCell',  {0},'ChildCell', {0} );

for i=1:m_new
    cordnew1(i).PixelList=cordnew(i).PixelList;
    cordnew1(i).PixelIdxList=cordnew(i).PixelIdxList;
    cordnew1(i).BoundingBox=cordnew(i).BoundingBox;
    cordnew1(i).Area=cordnew(i).Area;
    cordnew1(i).ConvexArea=cordnew(i).ConvexArea;
    cordnew1(i).EulerNumber=cordnew(i).EulerNumber;
    cordnew1(i).ConvexHull=cordnew(i).ConvexHull;
    cordnew1(i).Extrema=cordnew(i).Extrema;
    cordnew1(i).Centroid=cordnew(i).Centroid;
    cordnew1(i).Solidity=cordnew(i).Solidity;
    cordnew1(i).Extent=cordnew(i).Extent;
    cordnew1(i).EquivDiameter=cordnew(i).EquivDiameter;
    cordnew1(i).Perimeter=cordnew(i).Perimeter;
    cordnew1(i).MajorAxisLength=cordnew(i).MajorAxisLength;
    cordnew1(i).MinorAxisLength=cordnew(i).MinorAxisLength;
    cordnew1(i).Orientation=cordnew(i).Orientation;
    cordnew1(i).FilledArea=cordnew(i).FilledArea;
    cordnew1(i).SliceNo = 0;
    cordnew1(i).ParentCell=0;
    cordnew1(i).ChildCell=0;
end

for i=1:n
    c=int2str(i);
     text(cordnew(i).Centroid(1),cordnew(i).Centroid(2),c,'Color',[0,1,0],'FontSize',9);
     %plot(cordnew(i).Centroid(1),cordnew(i).Centroid(2),'Color',[1,0,0],'Marker','.','MarkerSize',25); hold on;
end

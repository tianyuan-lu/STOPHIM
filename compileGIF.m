for n=1:51
    mat = imread(sprintf('Cells_%d.png',n));
    im{n} = mat;
end
filename = 'K27MModel.gif'; % Specify the output file name
for idx = 1:51
    [A,map] = rgb2ind(im{idx},256);
    if idx == 1
        imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.1);
    else
        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.1);
    end
end
function out=imNormalize99(im)

im=double(im);
temp=sort(im(:),'descend');
th1=temp(round(length(im(:))/100));
th2=min(im(:));

out=(im-th2)/(th1-th2);
out(out>1)=1;



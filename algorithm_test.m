clc
clear

tic;

% single image simply test
% uint8StegPic = imread("img-test-set/lena.bmp");

% image-set test algorithm
%<test
keyPath = 'img-test-set/';
picsInfo = dir('img-test-set/*.*');
for i = 3: length(picsInfo)
    fprintf("%d. it is %s\n", i - 2, picsInfo(i, 1).name);
    uint8StegPic = imread(join([keyPath picsInfo(i, 1).name]));
%test/>

[picRow, picCol] = size(uint8StegPic);

%若列数为奇，则每列的最后一个像素不作隐藏用
picColUsed = picCol - mod(picCol, 2);

stegPic = double(uint8StegPic);

%在这里对隐写图片每一个像素都进行了加一，以对应MatLab中索引从1开始
stegPicSelfInc = stegPic + 1;

fiveBinSecrets = randi([0, 4], 1, picRow * picColUsed);

carrierPicDouble1 = stegPic;
carrierPicDouble2 = stegPic;

%构造魔法矩阵
global matSize;
matSize = 256;
global magicMat;
magicMat= zeros(matSize, matSize);
for i = 1: matSize
    for j = 1: matSize
        magicMat(i, j) = mod(j+1+(matSize-i-1)*2, 5);
    end
end

%这里要注意魔法矩阵不能使用的区域（像素值为0或255），注意判断
%另外要注意像素值已经加一，需在修改时减一
insertIndex = 1;
for i = 1: picRow
    for j = 1: 2: picColUsed
        %魔法矩阵不可使用的区域
        if isequal(stegPicSelfInc(i, j), 1) || isequal(stegPicSelfInc(i, j), 256)...
                || isequal(stegPicSelfInc(i, j+1), 1) || isequal(stegPicSelfInc(i, j+1), 256)
            continue;
        else
            %进行藏入
            
            curPixPair = [stegPicSelfInc(i, j) stegPicSelfInc(i, j+1)];
            %查找隐藏坐标
            
            %由于并不是每个像素对都可以插入密字，所以密字序号与i、j无关，所以需要单独设置
            src = crossSearch(curPixPair, fiveBinSecrets(1, insertIndex));
            dest = forkSearch(curPixPair, fiveBinSecrets(1, insertIndex+1));
            %歧义消除
            srcTrans = disambig(src, dest, curPixPair);
            %这里应该要减去一，并写入载体图片
            srcSelfDec = srcTrans - 1;
            destSelfDec = dest - 1;
            
            carrierPicDouble1(i, j) = srcSelfDec(1, 1);
            carrierPicDouble1(i, j+1) = srcSelfDec(1, 2);
            carrierPicDouble2(i, j) = destSelfDec(1, 1);
            carrierPicDouble2(i, j+1) = destSelfDec(1, 2);
            
            %密字序号更新
            insertIndex = insertIndex + 2;
        end
    end
end

insertSecretsNum = insertIndex - 1;

carrierPic1 = uint8(carrierPicDouble1);
carrierPic2 = uint8(carrierPicDouble2);

%提取
recPic1 = double(carrierPic1);
recPic2 = double(carrierPic2);

orginPicDouble = recPic1;
[picRow, picCol] = size(orginPicDouble);
picColUsed = picCol - mod(picCol, 2);

extractSecrets = zeros(1, picRow * picColUsed);
extractIndex = 1;

for i = 1: picRow
    for j = 1: 2: picColUsed
         if (isequal(recPic1(i, j), 0) || isequal(recPic1(i, j), 255) || isequal(recPic1(i, j+1), 0)...
                 || isequal(recPic1(i, j+1), 255)) && isequal([recPic1(i, j) recPic1(i, j+1)], [recPic2(i, j) recPic2(i, j+1)])
            continue;
         else
            src = [recPic1(i, j) recPic1(i, j+1)];
            dest = [recPic2(i, j) recPic2(i, j+1)];
            
            [center, isrc] = specRecog(src, dest);
            %恢复原图
            orginPicDouble(i, j)  = center(1, 1);
            orginPicDouble(i, j+1) = center(1, 2);
            
            %提取密字
            memIsrcSelfInc = math2mem(isrc + 1);
            memDestSelfInc = math2mem(dest + 1);
            extractSecrets(1, extractIndex) = magicMat(memIsrcSelfInc(1, 1), memIsrcSelfInc(1, 2));
            extractSecrets(1, extractIndex+1) = magicMat(memDestSelfInc(1, 1), memDestSelfInc(1, 2));
            
            extractIndex = extractIndex + 2;
         end
    end
end

extractSecretsNum = extractIndex - 1;

orginPic = uint8(orginPicDouble);

% 判断提取后的图片是否与原图一样
if isequal(uint8StegPic, orginPic)
    fprintf("pics are same\n");
else
    fprintf("pics are not same\n");
end

% 判断提取出的密文是否与藏入的密文一样
if  isequal(fiveBinSecrets(1: insertSecretsNum), extractSecrets(1: extractSecretsNum))
    fprintf("secret is same\n");
else
    fprintf("secret is not same\n");
end

capacity = (insertSecretsNum / (2 * picRow * picCol)) * log2(5);
fprintf("capacity is %f\n", capacity);

carrierPic1PSNR = psnr(carrierPic1, uint8StegPic);
carrierPic2PSNR = psnr(carrierPic2, uint8StegPic);
fprintf("carrierPic1PSNR is %f\n", carrierPic1PSNR);
fprintf("carrierPic2PSNR is %f\n", carrierPic2PSNR);

%<test
toc;
fprintf('\n');
end
%test/>

outputPics = {uint8StegPic carrierPic1 carrierPic2 orginPic};
for i = 1: 4
    subplot(2, 2, i);
    imshow(outputPics{1, i}, [256 256]);
end

toc;

%+形区域查找
function src = crossSearch(pixPair, secret)
    global matSize;
    global magicMat;
    crossDirect = [0 0; 0 1; 1 0; 0 -1; -1 0];
    [searchTime, ~] = size(crossDirect);
    memCoord = math2mem(pixPair);
    memCoordX = memCoord(1, 1);
    memCoordY = memCoord(1, 2);
    for c = 1: searchTime
        traversalNum = magicMat(memCoordX + crossDirect(c, 1), memCoordY + crossDirect(c, 2));
        if isequal(traversalNum, secret)
            mathCoord = mem2math([memCoordX + crossDirect(c, 1) memCoordY + crossDirect(c, 2)]);
            src = mathCoord;
            break;
        end
    end
end

%x形区域查找
function dest = forkSearch(pixPair, secret)
    global matSize;
    global magicMat;
    forkDirect = [0 0; -1 1; 1 1; 1 -1; -1 -1];
    [searchTime, ~] = size(forkDirect);
    memCoord = math2mem(pixPair);
    memCoordX = memCoord(1, 1);
    memCoordY = memCoord(1, 2);
    for f = 1: searchTime
        traversalNum = magicMat(memCoordX + forkDirect(f, 1), memCoordY + forkDirect(f, 2));
        if isequal(traversalNum, secret)
            mathCoord = mem2math([memCoordX + forkDirect(f, 1) memCoordY + forkDirect(f, 2)]);
            dest = mathCoord;
            break;
        end
    end
end

%歧义消除
function srcTrans = disambig(src, dest, center)
    srcTrans = src;
    s2d = dest - src;
    lenSqu = s2d(1, 1).^2 + s2d(1, 2).^2;
    c2s = src - center;
    c2d = dest - center;
    if isequal(lenSqu, 1)
        forkMultiRes = forkMulti(c2s, c2d);
        %顺时针
        if forkMultiRes > 0
            srcTrans = src - s2d;
        %逆时针
        elseif forkMultiRes < 0
            srcTrans = dest - 2 * c2d;
        end
    end
end

%判断vec2与vec1的位置关系
function forkMultiRes = forkMulti(vec1, vec2)
    %x1*y2 - x2*y1
    forkMultiRes =  vec2(1, 1) * vec1(1, 2) - vec2(1, 2) * vec1(1, 1);
end

%种类识别
function [center, isrc] = specRecog(src, dest)
    %矢量长度平方为0
    if isequal(src, dest)
        center = src;
        isrc = src;
    else
        lenSqu = (dest(1, 1) - src(1, 1)).^2 + (dest(1, 2) - src(1, 2)).^2;
        %矢量长度平方为1
        if isequal(lenSqu, 1)
            center = dest;
            isrc = src;
        %矢量长度平方为2
        elseif isequal(lenSqu, 2)
            center = src;
            isrc = src;
        %矢量长度平方为5
        elseif isequal(lenSqu, 5)
            vec = [dest(1, 1) - src(1, 1) dest(1, 2)-src(1, 2)];
            center = src + fix(vec / 2);
            isrc = src;
        %矢量长度平方为8
        elseif isequal(lenSqu, 8)
             center = [floor((src(1, 1) + dest(1, 1)) / 2) floor((src(1, 2) + dest(1, 2)) / 2)];
             vec = [dest(1, 1) - src(1, 1) dest(1, 2) - src(1, 2)];
             if isequal(vec, [-2, 2])
                 isrc = dest + [1 0];
             elseif isequal(vec, [2 2])
                 isrc = dest + [0 -1];
             elseif isequal(vec, [2 -2])
                 isrc = dest + [-1 0];
             elseif isequal(vec, [-2 -2])
                 isrc = dest + [0 1];
             end
        %矢量长度平方为4
        elseif isequal(lenSqu, 4)
            vec = [dest(1, 1) - src(1, 1) dest(1, 2) - src(1, 2)];
            if isequal(vec, [2 0])
                center = src + [1 -1];
            elseif isequal(vec, [0 -2])
                center = src + [-1 -1];
            elseif isequal(vec, [-2 0])
                center = src + [-1 1];
            elseif isequal(vec, [0 2])
                center = src + [1 1];
            end
            isrc = src + [vec(1, 1) / 2 vec(1, 2) / 2];
        end
    end
end

function memCoord = math2mem(mathCoord)
    %(i, j)->(matSize+1- j, i)
    global matSize;
    memCoord = [matSize + 1 - mathCoord(1, 2) mathCoord(1, 1)];
end

function mathCoord = mem2math(memCoord)
    %(i, j) -> (j, matSize + 1 - i)
    global matSize;
    mathCoord = [memCoord(1, 2) matSize + 1 - memCoord(1, 1)];
end
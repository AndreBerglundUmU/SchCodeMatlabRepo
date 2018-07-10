function appendToText(fileName,M,N,sigma,P,PLow,PUpp)
    fileID = fopen(fileName,'a');
    fprintf(fileID,'%5d %5d %5d %5d %5d %5d\n',[M N sigma P PLow PUpp]);
    fclose(fileID);
end
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main() {

// I/O Managment 
    FILE *port = fopen("COM1", "w+"); // Replace "COM1" with the actual port name
    
    if (port == NULL) {
        printf("Error opening port!\n");
        return 1;
    }

// Memory Managment
    char data_in0[1000];
    char *data_in1 = (char *)malloc(1000 * sizeof(char)); 

    fflush(port); // Flush the buffer to ensure data is sent immediately
    fread(&data_in1+1, sizeof(char), 1, port);

    
// data transfer
memcpy(data_in1, data_in0, 1000);


// Data Processing
for (int i=0; i<990 ; i++)
data_in1[i]=data_in1[i]+data_in1[i+1]+data_in1[i+2]+data_in1[i+3]+data_in1[i+4]+data_in1[i+5]+data_in1[i+6]+data_in1[i+7]+data_in1[i+8]+data_in1[i+9];
    printf("Received Data: %s\n", data_in1);


memcpy(data_in0, data_in1, 1000);
    fclose(port);

    return 0;
}

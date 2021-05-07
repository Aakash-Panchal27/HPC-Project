/*
Commands to compile and run this code:

Compile: g++ -o test.o -fopenmp -std=gnu++11 serial.cpp
Run: ./test.o
*/

#include<bits/stdc++.h>
#include<fstream>
using namespace std;

typedef std::basic_string<unsigned char> ustring;

const int bits_unit = 8;

struct node {
    int frequency;
    unsigned char in; // 4 byte or 8 byte
    node * left;
    node * right;
    node(int freq, unsigned char i)
    {
        frequency = freq;
        in = i;
        left = NULL;
        right = NULL;
    }
   
};

struct comparator{

    bool operator()(node* temp1, node* temp2)
    {
        if(temp1->frequency > temp2->frequency)
            return true;
        else 
            return false;
    }
};

void obtain_huffmancodes(node* root, int byte, pair<int, int>* encode, int cur_bit_id = 0)
{
    if(root->left == NULL && root->right == NULL)
    {
        encode[root->in].first = byte;
        encode[root->in].second = cur_bit_id;
        assert(cur_bit_id <= 32);
        return;
    }
    int nxt_id = cur_bit_id + 1;
    int nxt_byte = byte | (1 << cur_bit_id);
    obtain_huffmancodes(root->left, byte, encode, nxt_id);
    obtain_huffmancodes(root->right, nxt_byte, encode, nxt_id);
}

int main()
{
    cout << "Enter the file name with extension (e.g. bonsai.raw)" << endl;
    string filename;
    cin >> filename;
    string decompressed_file_name = "output_" + filename;

    // Stores the no. of times we will repeat the whole procedure to measure the time accurately by taking average in the end.
    long double runs = 10;
    // To store times for different steps
    long double freq_table_time,huffman_tree_time,hashtable_time,compression_time,decompression_time,whole_time;
    freq_table_time = huffman_tree_time = hashtable_time = compression_time = decompression_time = whole_time = 0;
    bool printed = false;

    for(int run = 0; run < runs; run++) {

        struct timespec whole_start, whole_end;
        clock_gettime(CLOCK_MONOTONIC, &whole_start);
        
        // Read file into buffer and find its size
        ifstream file(filename, ios::in | ios::binary);
        file.seekg (0, file.end);
        int file_size = file.tellg();
        assert(file_size <= (1 << 29));
        file.seekg (0, file.beg);
        unsigned char* buffer = new unsigned char[file_size];
        file.read((char*)buffer, file_size);
        file.close();
        // Create frequency table--------------------------------------------------
        auto start = chrono::high_resolution_clock::now();

        int freq[256] = {0};
        for(int i = 0; i < file_size; i++)
            freq[buffer[i]]++;
            
        auto end = chrono::high_resolution_clock::now();
        double total_time =  
          chrono::duration_cast<chrono::nanoseconds>(end - start).count();
        total_time *= 1e-9;
        freq_table_time += total_time;
        // Create frequency table--------------------------------------------------

        // Create huffman tree-----------------------------------------------------
        start = chrono::high_resolution_clock::now();

        priority_queue<node*, vector<node*>, comparator> container;

        // Standard algo
        for(int i = 0; i < 256; i++)
        {
            if(freq[i] > 0) {
                node* temp = new node(freq[i], i);
                container.push(temp);
            }
        }

        while(container.size() > 1)
        {
            node* left = container.top();
            container.pop();
            
            node * right = container.top();
            container.pop();
            
            node* curr = new node(left->frequency + right->frequency, '0');
            curr->left = left;
            curr->right = right;
            container.push(curr);
        }

        node* root = container.top();
        container.pop();
        node* root_store = root;

        end = chrono::high_resolution_clock::now();
        total_time =  
          chrono::duration_cast<chrono::nanoseconds>(end - start).count();
        total_time *= 1e-9;
        huffman_tree_time += total_time;
        // Create huffman tree-------------------------------------------------

        // Creating Hashtable of encodes----------------------------------------
        start = chrono::high_resolution_clock::now();

        // Create hashtable so that we don't have to traverse the tree again and again
        int ans = 0;
        pair<int, int> encode[256];
        obtain_huffmancodes(root, ans, encode);
        
        end = chrono::high_resolution_clock::now();
        total_time =  
          chrono::duration_cast<chrono::nanoseconds>(end - start).count();
        total_time *= 1e-9;
        hashtable_time += total_time;
        // Creating Hashtable of encodes----------------------------------------

        // Compress--------------------------------------------------
        start = chrono::high_resolution_clock::now();

        int idx = 7;
        ustring comp_file;
        unsigned int comp_total_bits = 0;
        unsigned char a;
        unsigned char c = 0;
        for(int i = 0; i < file_size; i++) {
            a = buffer[i];
            int byte = encode[a].first;
            int len = encode[a].second;
            comp_total_bits += len;
            for(int i = 0; i < len; i++){
                bool bit = byte & (1 << i);
                c |= (bit << idx);
                idx--;
                if(idx < 0) {
                    comp_file += c;
                    idx = 7;
                    // refresh c
                    c &= 0;
                }
            }
        }

        if(idx != 7)
            comp_file += c;
        
        if(!printed) {
            cout << "Compressed file size in bytes: " << comp_file.size() << endl;
            printed = true;
        }

        end = chrono::high_resolution_clock::now();
        total_time =  
          chrono::duration_cast<chrono::nanoseconds>(end - start).count();
        total_time *= 1e-9;
        compression_time += total_time;

        delete[] buffer;
        // Compress----------------------------------------------------
        
        // Decompress---------------------------------------------------
        start = chrono::high_resolution_clock::now();
        
        int comp_file_size = comp_file.length();
        unsigned int total_bits = comp_total_bits;
        
        ofstream decompressed(decompressed_file_name, ios::out | ios::binary);
        buffer = new unsigned char[file_size];
        int byte_no = 0;
        root = root_store;
        for(int i = 0; i < comp_file_size && total_bits > 0; i++) {
            a = comp_file[i];
            idx = 7;
            while(idx >= 0 && total_bits > 0) {
                if(a & (1 << idx))
                    root = root->right;
                else
                    root = root->left;
                if(!root->right && !root->left) {
                    buffer[byte_no++] = root->in;
                    root = root_store;
                }
                idx--;
                total_bits--;
            }
        }
        
        end = chrono::high_resolution_clock::now();
        total_time =  
          chrono::duration_cast<chrono::nanoseconds>(end - start).count();
        total_time *= 1e-9;
        decompression_time += total_time;
        
        decompressed.write((char*)buffer, byte_no);
        decompressed.close();
        delete[] buffer;
        // Decompress-----------------------------------------------------------

        clock_gettime(CLOCK_MONOTONIC, &whole_end);
        total_time = (whole_end.tv_sec - whole_start.tv_sec) * 1e9;
        total_time = (total_time + (whole_end.tv_nsec - whole_start.tv_nsec)) * 1e-9;
        whole_time += total_time;
    }

    freq_table_time /= double(runs);
    huffman_tree_time /= double(runs);
    hashtable_time /= double(runs);
    compression_time /= double(runs);
    decompression_time /= double(runs);
    whole_time /= double(runs);

    cout << fixed << setprecision(20) << "Creating frequency table: " << freq_table_time << endl;
    cout << fixed << setprecision(20) << "Creating Huffman Tree: " << huffman_tree_time << endl;
    cout << fixed << setprecision(20) << "Creating hashtable of encodes: " << hashtable_time << endl;
    cout << fixed << setprecision(20) << "Compression: " << compression_time << endl;
    cout << fixed << setprecision(20) << "Decompression: " << decompression_time << endl;
    cout << fixed << setprecision(20) << "Total time: " << whole_time << endl;
    return 0;
}
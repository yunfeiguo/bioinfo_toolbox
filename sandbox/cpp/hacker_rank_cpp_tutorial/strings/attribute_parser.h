//
// Created by guoy28 on 10/24/17.
//

#ifndef HACKER_RANK_CPP_TUTORIAL_ATTRIBUTE_PARSER_H
#define HACKER_RANK_CPP_TUTORIAL_ATTRIBUTE_PARSER_H

#endif //HACKER_RANK_CPP_TUTORIAL_ATTRIBUTE_PARSER_H
#include <iostream>
#include <cmath>
#include <cstdio>
#include <vector>
#include <algorithm>
#include <stack>
#include <map>
#include <set>
using namespace std;

class AttributeParser {
public:
    AttributeParser(std::istream*, std::ostream*);
};

AttributeParser::AttributeParser(std::istream* in, std::ostream* out) {
    int N,Q;
    *in >> N >> Q;
    stack<std::string> current_tags;
    map<string, map<string, string> > tag2attribute_map;
    map<string, set<string> > parent2child; //parent tag - child tag
    map<string, set<string> > child2parent; //child tag - parent tag
    for (int i = 0; i < N; i++) {
//read input
        std::string token;
        *in >> token;
        if (token[1] == '/') { //this is a closing tag
//assume: current_tags.top() == tag.substr(2,tag.size()-3);
            current_tags.pop();
        } else {
            string tag = token.substr(1, token.size()- (token[token.size()-1] == '>' ? 2 : 1));
            if (!current_tags.empty()) {
                parent2child[current_tags.top()].insert(tag);
                child2parent[tag].insert(current_tags.top());
            }
            current_tags.push(tag); //assume each starting tag only occurs once
            if (tag2attribute_map.count(current_tags.top()) == 0) {
                map<string, string> attribute2value;
                tag2attribute_map[current_tags.top()] = attribute2value;
            }

            while(token[token.size() - 1] != '>') {
                *in >> token;
                string attribute = token;
                *in >> token; //=
                *in >> token;
                string value = token.substr(1, token.size() - (token[token.size() - 1] == '>' ? 3 : 2));
                tag2attribute_map[current_tags.top()][attribute] = value;
            }
        }
    }
    for (int i = 0; i < Q; i++) {
        string token;
        string previous_tag;
        string current_tag;
        *in >> token;
        int previous_dot_position = -1;
        size_t dot_position = token.find('.');
        size_t tilde_position = token.find_last_of('~');
        bool not_found = false;
        while ((previous_dot_position == -1 && dot_position == string::npos) || previous_dot_position != dot_position) {
            current_tag = token.substr(previous_dot_position + 1, (dot_position == string::npos? tilde_position : dot_position) - previous_dot_position - 1);
            if (!previous_tag.empty()) {
                if (parent2child.count(previous_tag) == 0 || parent2child[previous_tag].count(current_tag) == 0) {
                    not_found = true;
                    break;
                }
            } else {
                if (child2parent.count(current_tag) != 0) {
                    not_found = true;
                    break;
                }
            }
            previous_tag = current_tag;
            previous_dot_position = dot_position;
            if (dot_position == string::npos) {
                break;
            }
            dot_position = token.find('.', dot_position + 1);
        }
        if (not_found) {
            *out << "Not Found!" << endl;
            continue;
        }
        string tag;
        dot_position = token.find_last_of('.');
        if (dot_position == string::npos) {
            tag = token.substr(0, tilde_position);
        } else {
            tag = token.substr(dot_position + 1, tilde_position - dot_position - 1);
        }
        string attribute = token.substr(tilde_position + 1, token.size() - 1);
        if (tag2attribute_map.count(tag) == 0 || tag2attribute_map[tag].count(attribute) == 0) {
            *out << "Not Found!" << endl;
            continue;
        }
        *out << tag2attribute_map[tag][attribute] << endl;
    }
}

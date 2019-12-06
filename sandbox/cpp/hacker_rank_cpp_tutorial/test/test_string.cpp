//
// Created by guoy28 on 10/24/17.
//

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_TOOLS_UNDER_DEBUGGER
#define BOOST_TEST_MODULE test_string
#include <boost/test/unit_test.hpp>
#include <string>
#include <boost/filesystem.hpp>
#include <fstream>
#include <iostream>
#include "../strings/attribute_parser.h"

BOOST_AUTO_TEST_CASE(simple_io) {
    boost::filesystem::path temp = boost::filesystem::unique_path();
    const std::string temp_name = temp.native();
    std::ifstream ifs("../test/test_files/attribute_parser_example_input.txt", std::ifstream::in);
    boost::filesystem::ofstream ofs(temp, std::ofstream::out);
    AttributeParser a(&ifs, &ofs);
    ofs.close();

    boost::filesystem::ifstream test_output(temp, std::ofstream::in);
    std::ifstream expected_output("../test/test_files/attribute_parser_example_output.txt", std::ofstream::in);

    std::istream_iterator<char> test_output_itr_begin(test_output), test_output_itr_end;
    std::istream_iterator<char> expected_output_itr_begin(expected_output), expected_output_itr_end;

    BOOST_CHECK_EQUAL_COLLECTIONS(test_output_itr_begin, test_output_itr_end, expected_output_itr_begin, expected_output_itr_end);
}

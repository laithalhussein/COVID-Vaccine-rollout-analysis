function [US_full_ID, US_abb_ID] = extract_suppl_data

US_full_ID = {
'Alabama'
'Alaska'
'American Samoa'
'Arizona'
'Arkansas'
'California'
'Canal Zone'
'Colorado'
'Connecticut'
'Delaware'
'District of Columbia'
'Florida'
'Georgia'
'Guam'
'Hawaii'
'Idaho'
'Illinois'
'Indiana'
'Iowa'
'Kansas'
'Kentucky'
'Louisiana'
'Maine'
'Mariana Islands'
'Marshall Islands'
'Maryland'
'Massachusetts'
'Michigan'
'Micronesia' 
'Minnesota'
'Mississippi'
'Missouri'
'Montana'
'Nebraska'
'Nevada'
'New Hampshire'
'New Jersey'
'New Mexico'
'New York'
'North Carolina'
'North Dakota'
'Ohio'
'Oklahoma'
'Oregon'
'Palau'
'Pennsylvania'
'Puerto Rico'
'Rhode Island'
'South Carolina'
'South Dakota'
'Tennessee'
'Texas'
'Utah'
'Vermont'
'U.S. Virgin Islands'
'Virginia'
'Washington'
'West Virginia'
'Wisconsin'
'Wyoming'
'Federal Entities'
};

US_abb_ID = {
'AL'
'AK'
'AS'
'AZ'
'AR'
'CA'
'CZ'
'CO'
'CT'
'DE'
'DC'
'FL'
'GA'
'GU'
'HI'
'ID'
'IL'
'IN'
'IA'
'KS'
'KY'
'LA'
'ME'
'MP'
'MH'
'MD'
'MA'
'MI'
'FM'
'MN'
'MS'
'MO'
'MT'
'NE'
'NV'
'NH'
'NJ'
'NM'
'NY'
'NC'
'ND'
'OH'
'OK'
'OR'
'RP'
'PA'
'PR'
'RI'
'SC'
'SD'
'TN'
'TX'
'UT'
'VT'
'VI'
'VA'
'WA'
'WV'
'WI'
'WY'
'US'
};

%save('US_ID_names.mat', 'US_full_ID', 'US_abb_ID');
return
